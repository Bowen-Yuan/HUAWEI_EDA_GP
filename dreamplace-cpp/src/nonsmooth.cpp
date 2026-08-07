#include "nonsmooth.h"

#include "bookshelf.h"
#include "wirelength.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace dpcpp {
namespace {

struct BundleCut {
    Real hpwl = 0.0;
    std::vector<Real> x, y, gx, gy;
};

struct BundleStats {
    int cuts = 0;
    Real newest_weight = 1.0;
    Real model_error = 0.0;
    Real cosine = 1.0;
};

std::vector<Real> capture(const Database& db, const std::vector<Filler>& fillers) {
    const std::size_t n = db.movable_ids.size() + fillers.size();
    std::vector<Real> values(2 * n);
    std::size_t k = 0;
    for (int id : db.movable_ids) values[k++] = db.nodes[id].x;
    for (const Filler& filler : fillers) values[k++] = filler.x;
    k = n;
    for (int id : db.movable_ids) values[k++] = db.nodes[id].y;
    for (const Filler& filler : fillers) values[k++] = filler.y;
    return values;
}

void apply(Database& db, std::vector<Filler>& fillers,
           const std::vector<Real>& values) {
    const std::size_t n = db.movable_ids.size() + fillers.size();
    if (values.size() != 2 * n) throw std::runtime_error("position vector size mismatch");
    std::size_t k = 0;
    for (int id : db.movable_ids) db.nodes[id].x = values[k++];
    for (Filler& filler : fillers) filler.x = values[k++];
    k = n;
    for (int id : db.movable_ids) db.nodes[id].y = values[k++];
    for (Filler& filler : fillers) filler.y = values[k++];
    clamp_to_region(db, fillers);
}

Real l1_norm(const std::vector<Real>& values) {
    Real sum = 0.0;
    for (Real value : values) sum += std::abs(value);
    return sum;
}

Real rms_movable(const Database& db, const std::vector<Real>& gx,
                 const std::vector<Real>& gy) {
    Real sum = 0.0;
    for (int id : db.movable_ids) sum += gx[id] * gx[id] + gy[id] * gy[id];
    return std::sqrt(sum / std::max<std::size_t>(1, 2 * db.movable_ids.size()));
}

void add_gp_noise(Database& db, std::vector<Filler>& fillers,
                  std::uint64_t seed, Real ratio) {
    if (ratio <= 0.0) return;
    std::mt19937_64 generator(seed ^ 0xd1b54a32d192ed03ULL);
    std::uniform_real_distribution<Real> unit(-0.5, 0.5);
    for (int id : db.movable_ids) {
        Node& node = db.nodes[id];
        node.x += unit(generator) * node.width * ratio;
        node.y += unit(generator) * node.height * ratio;
    }
    for (Filler& filler : fillers) {
        filler.x += unit(generator) * filler.width * ratio;
        filler.y += unit(generator) * filler.height * ratio;
    }
    clamp_to_region(db, fillers);
}

Real vector_cosine(const Database& db, const std::vector<Real>& ax,
                   const std::vector<Real>& ay, const BundleCut& b) {
    Real dot = 0.0, an = 0.0, bn = 0.0;
    for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
        const int id = db.movable_ids[i];
        dot += ax[id] * b.gx[i] + ay[id] * b.gy[i];
        an += ax[id] * ax[id] + ay[id] * ay[id];
        bn += b.gx[i] * b.gx[i] + b.gy[i] * b.gy[i];
    }
    return dot / std::max<Real>(std::sqrt(an * bn), 1.0e-30);
}

BundleStats apply_bundle(const Database& db, Real hpwl, Real learning_rate,
                         int iteration, const GlobalPlaceConfig& config,
                         std::vector<BundleCut>& cuts, int& last_cut,
                         std::vector<Real>& gx, std::vector<Real>& gy) {
    BundleStats stats;
    Real cosine = cuts.empty() ? 1.0 : vector_cosine(db, gx, gy, cuts.back());
    const int age = iteration - last_cut;
    bool insert = cuts.empty() || age >= config.bundle_interval;
    if (!cuts.empty() && age >= config.bundle_min_interval &&
        cosine <= config.bundle_cosine_trigger) insert = true;
    if (insert) {
        BundleCut cut;
        cut.hpwl = hpwl;
        const std::size_t n = db.movable_ids.size();
        cut.x.resize(n); cut.y.resize(n); cut.gx.resize(n); cut.gy.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            const int id = db.movable_ids[i];
            cut.x[i] = db.nodes[id].x; cut.y[i] = db.nodes[id].y;
            cut.gx[i] = gx[id]; cut.gy[i] = gy[id];
        }
        if (static_cast<int>(cuts.size()) >= config.bundle_size) cuts.erase(cuts.begin());
        cuts.push_back(std::move(cut));
        last_cut = iteration;
    }

    const int kcount = static_cast<int>(cuts.size());
    stats.cuts = kcount;
    stats.cosine = cosine;
    if (kcount == 0) return stats;
    std::vector<Real> beta(kcount, 0.0), gram(kcount * kcount, 0.0);
    for (int j = 0; j < kcount; ++j) {
        Real shift = 0.0;
        for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
            const int id = db.movable_ids[i];
            shift += cuts[j].gx[i] * (db.nodes[id].x - cuts[j].x[i]) +
                     cuts[j].gy[i] * (db.nodes[id].y - cuts[j].y[i]);
        }
        beta[j] = cuts[j].hpwl + shift;
        for (int k = 0; k <= j; ++k) {
            Real value = 0.0;
            for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
                value += cuts[j].gx[i] * cuts[k].gx[i] +
                         cuts[j].gy[i] * cuts[k].gy[i];
            }
            gram[j * kcount + k] = gram[k * kcount + j] = value;
        }
    }
    const Real tau = config.bundle_prox_scale * learning_rate /
                     std::max<Real>(rms_movable(db, gx, gy), 1.0e-9);
    std::vector<Real> alpha(kcount, 0.0);
    alpha.back() = 1.0;
    for (int fw = 0; fw < 20; ++fw) {
        std::vector<Real> dual_gradient(kcount, 0.0);
        int best = 0;
        for (int j = 0; j < kcount; ++j) {
            Real q = 0.0;
            for (int k = 0; k < kcount; ++k) q += gram[j * kcount + k] * alpha[k];
            dual_gradient[j] = beta[j] - tau * q;
            if (dual_gradient[j] > dual_gradient[best]) best = j;
        }
        std::vector<Real> direction(kcount);
        for (int j = 0; j < kcount; ++j) direction[j] = -alpha[j];
        direction[best] += 1.0;
        Real numerator = 0.0, curvature = 0.0;
        for (int j = 0; j < kcount; ++j) {
            numerator += dual_gradient[j] * direction[j];
            for (int k = 0; k < kcount; ++k)
                curvature += direction[j] * gram[j * kcount + k] * direction[k];
        }
        if (numerator <= 1.0e-7) break;
        const Real fraction = curvature > 1.0e-20
            ? std::min<Real>(1.0, numerator / (tau * curvature)) : 1.0;
        for (int j = 0; j < kcount; ++j) alpha[j] += fraction * direction[j];
    }

    std::vector<Real> aggregate_x(db.nodes.size(), 0.0);
    std::vector<Real> aggregate_y(db.nodes.size(), 0.0);
    Real aggregate_beta = 0.0;
    for (int j = 0; j < kcount; ++j) {
        aggregate_beta += alpha[j] * beta[j];
        for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
            const int id = db.movable_ids[i];
            aggregate_x[id] += alpha[j] * cuts[j].gx[i];
            aggregate_y[id] += alpha[j] * cuts[j].gy[i];
        }
    }
    const Real mix = std::clamp(config.bundle_current_mix, 0.0, 1.0);
    for (int id : db.movable_ids) {
        gx[id] = mix * gx[id] + (1.0 - mix) * aggregate_x[id];
        gy[id] = mix * gy[id] + (1.0 - mix) * aggregate_y[id];
    }
    stats.newest_weight = alpha.back();
    stats.model_error = std::max<Real>(0.0, hpwl - aggregate_beta);
    return stats;
}

Real update_lambda_control(Real eta, LambdaPolicy policy, int density_step,
                           Real hpwl, Real overflow, Real previous_hpwl,
                           Real initial_overflow, Real& previous_overflow,
                           Real& integral_error,
                           const GlobalPlaceConfig& config) {
    if (policy == LambdaPolicy::Dreamplace) {
        const Real delta = hpwl - previous_hpwl;
        Real multiplier = 1.0;
        if (delta < 0.0) {
            multiplier = 1.05 * std::max(std::pow(0.9999, density_step), 0.98);
        } else {
            multiplier = std::clamp(
                std::pow(1.05, 1.0 - delta / config.reference_hpwl_delta),
                0.95, 1.05);
        }
        return std::clamp(eta * multiplier,
                          config.lambda_control_min, config.lambda_control_max);
    }
    if (config.lambda_update_interval <= 0 ||
        density_step % config.lambda_update_interval != 0) return eta;
    if (policy == LambdaPolicy::Ratio) {
        const Real hpwl_target = config.hpwl_baseline > 0.0
            ? config.hpwl_baseline : hpwl;
        const Real hpwl_ratio = hpwl / std::max<Real>(hpwl_target, 1.0);
        const Real overflow_ratio = overflow /
                                    std::max<Real>(config.overflow_baseline, 1.0e-9);
        Real factor = 1.0;
        // Lambda multiplies density: excess overflow must increase it, while
        // an HPWL-dominated feasible state can release density pressure.
        if (overflow_ratio > 1.05 * hpwl_ratio) factor = 1.5;
        else if (hpwl_ratio > 1.05 * overflow_ratio) factor = 2.0 / 3.0;
        if (overflow > 1.05 * config.stop_overflow) factor = std::max(factor, 1.25);
        return std::clamp(eta * factor,
                          config.lambda_control_min, config.lambda_control_max);
    }

    const Real progress = std::min<Real>(1.0, density_step /
        static_cast<Real>(std::max(1, config.lambda_trajectory_horizon)));
    const Real smoothstep = progress * progress * (3.0 - 2.0 * progress);
    const Real desired = config.stop_overflow +
        (initial_overflow - config.stop_overflow) * (1.0 - smoothstep);
    const Real scale = std::max<Real>(config.stop_overflow, 1.0e-4);
    const Real error = std::clamp((overflow - desired) / scale, -2.0, 4.0);
    const Real desired_drop = (previous_overflow - desired) / scale;
    const Real actual_drop = (previous_overflow - overflow) / scale;
    const Real trend_error = std::clamp(desired_drop - actual_drop, -2.0, 2.0);
    if (progress < 1.0) integral_error = std::clamp(integral_error + error, -3.0, 3.0);
    else integral_error = 0.0;
    // The feedback is evaluated only once per interval.  Express the gains
    // as an equivalent per-iteration log change and accumulate them across
    // the interval; otherwise a 25-step controller ramps about 25x too
    // slowly compared with DREAMPlace's per-iteration lambda update.
    Real log_delta_per_iteration = 0.045 + 0.08 * error +
        0.003 * integral_error + 0.04 * trend_error;
    if (overflow <= config.stop_overflow && config.hpwl_baseline > 0.0) {
        const Real pressure = std::max<Real>(0.0, hpwl / config.hpwl_baseline - 1.0);
        log_delta_per_iteration -= 0.01 * pressure;
    }
    log_delta_per_iteration = std::clamp(
        log_delta_per_iteration, -0.05, 0.06);
    const Real log_delta = std::clamp(
        log_delta_per_iteration * config.lambda_update_interval, -1.0, 1.25);
    previous_overflow = overflow;
    return std::clamp(eta * std::exp(log_delta),
                      config.lambda_control_min, config.lambda_control_max);
}

Real coordinate_step(GlobalOptimizer optimizer, Real gradient, Real learning_rate,
                     int age, const GlobalPlaceConfig& config, Real& first,
                     Real& second, Real& maximum_second) {
    if (optimizer == GlobalOptimizer::HeavyBall) {
        first = config.momentum * first + (1.0 - config.momentum) * gradient;
        return learning_rate * first;
    }
    if (optimizer == GlobalOptimizer::AdaGrad) {
        second += gradient * gradient;
        return std::clamp(learning_rate * gradient /
            (std::sqrt(second) + config.epsilon),
            -config.max_step_multiplier * learning_rate,
             config.max_step_multiplier * learning_rate);
    }
    first = config.beta1 * first + (1.0 - config.beta1) * gradient;
    second = config.beta2 * second + (1.0 - config.beta2) * gradient * gradient;
    const Real first_hat = first / std::max<Real>(1.0e-12, 1.0 - std::pow(config.beta1, age));
    const Real second_hat = second / std::max<Real>(1.0e-12, 1.0 - std::pow(config.beta2, age));
    Real denominator = second_hat;
    if (optimizer == GlobalOptimizer::AMSGrad) {
        maximum_second = std::max(maximum_second, second_hat);
        denominator = maximum_second;
    }
    return std::clamp(learning_rate * first_hat /
        (std::sqrt(denominator) + config.epsilon),
        -config.max_step_multiplier * learning_rate,
         config.max_step_multiplier * learning_rate);
}

Real feasible_band_lambda(Real eta, Real overflow,
                          const GlobalPlaceConfig& config) {
    const Real upper = config.stop_overflow;
    const Real lower = config.refinement_lower_overflow;
    const Real band = std::max<Real>(upper - lower, 1.0e-4);
    Real log_delta = 0.0;
    if (overflow > upper) {
        log_delta = config.refinement_lambda_gain * (overflow - upper) / band;
    } else if (overflow < lower) {
        log_delta = -config.refinement_lambda_gain * (lower - overflow) / band;
    } else {
        // Release density pressure gradually inside the feasible band.  The
        // release vanishes at its upper edge, so an infeasible step reverses
        // the controller immediately instead of allowing lambda to run away.
        log_delta = -0.05 * config.refinement_lambda_gain *
                    (upper - overflow) / band;
    }
    log_delta = std::clamp(log_delta, -0.05, 0.05);
    return std::clamp(eta * std::exp(log_delta),
                      config.lambda_control_min, config.lambda_control_max);
}

void save_global_snapshot(const Database& db, const GlobalPlaceConfig& config,
                          int iteration) {
    if (config.snapshot_every <= 0 || config.snapshot_dir.empty() ||
        iteration % config.snapshot_every != 0) return;
    std::ostringstream name;
    name << "global_" << std::setw(4) << std::setfill('0') << iteration << ".pl";
    write_bookshelf_pl(db, (std::filesystem::path(config.snapshot_dir) /
                            name.str()).string());
}

}  // namespace

GlobalPlaceResult nonsmooth_global_place(Database& db,
                                         std::vector<Filler>& fillers,
                                         const GlobalPlaceConfig& config) {
    const auto start = std::chrono::steady_clock::now();
    ElectricDensity density(db, config.bins_x, config.bins_y, config.target_density);
    add_gp_noise(db, fillers, config.seed, config.gp_noise_ratio);
    std::ofstream metrics;
    if (!config.metrics_path.empty()) {
        std::filesystem::create_directories(
            std::filesystem::path(config.metrics_path).parent_path());
        metrics.open(config.metrics_path);
        metrics << "iteration,exact_hpwl,overflow,max_density,density_energy,"
                   "lambda_base,lambda_control,lambda_effective,learning_rate,"
                   "bundle_cuts,bundle_newest_weight,bundle_model_error,bundle_cosine,"
                   "refinement,optimizer\n";
        metrics << std::setprecision(12);
    }

    const std::size_t n = db.movable_ids.size() + fillers.size();
    std::vector<Real> first(2 * n, 0.0), second(2 * n, 0.0), maximum_second(2 * n, 0.0);
    std::vector<BundleCut> cuts;
    int last_cut = -1000000;
    Real lambda_base = 0.0;
    Real lambda_control = 1.0;
    Real previous_hpwl = exact_hpwl(db);
    Real initial_density_overflow = -1.0;
    Real previous_overflow = 0.0;
    Real integral_error = 0.0;
    int density_step = 0;
    bool refinement_active = false;
    int refinement_start = -1;
    GlobalPlaceResult result;
    std::vector<Real> best_feasible, best_overflow = capture(db, fillers);
    Metrics best_overflow_metrics;
    best_overflow_metrics.overflow = std::numeric_limits<Real>::infinity();

    for (int iteration = 0; iteration < config.iterations; ++iteration) {
        std::vector<Real> dgx, dgy, fgx, fgy;
        const DensityResult d = density.compute(db, fillers, &dgx, &dgy, &fgx, &fgy);
        std::vector<Real> wgx, wgy;
        const Real hpwl = exact_hpwl_subgradient(db, config.degree_limit, &wgx, &wgy);
        ++result.objective_evaluations;
        const bool density_active = iteration >= config.hpwl_only_iterations;
        if (density_active && lambda_base == 0.0) {
            Real wire_l1 = 0.0, density_l1 = 0.0;
            for (int id : db.movable_ids) {
                wire_l1 += std::abs(wgx[id]) + std::abs(wgy[id]);
                density_l1 += std::abs(dgx[id]) + std::abs(dgy[id]);
            }
            density_l1 += l1_norm(fgx) + l1_norm(fgy);
            lambda_base = config.density_weight_scale * wire_l1 /
                          std::max<Real>(density_l1, 1.0e-30);
            initial_density_overflow = d.overflow;
            previous_overflow = d.overflow;
            std::cout << "[Exact] lambda_base=" << lambda_base
                      << " wire_l1=" << wire_l1
                      << " density_l1=" << density_l1 << '\n';
        }
        if (config.feasible_refinement && density_active && !refinement_active &&
            d.overflow <= config.stop_overflow) {
            refinement_active = true;
            refinement_start = iteration;
            std::fill(first.begin(), first.end(), 0.0);
            std::fill(second.begin(), second.end(), 0.0);
            std::fill(maximum_second.begin(), maximum_second.end(), 0.0);
            cuts.clear();
            last_cut = -1000000;
            integral_error = 0.0;
            std::cout << "[Refine] start iter=" << iteration
                      << " hpwl=" << hpwl << " overflow=" << d.overflow
                      << " optimizer=" << global_optimizer_name(config.refinement_optimizer)
                      << " lr_scale=" << config.refinement_learning_rate_scale << '\n';
        }
        const Real lambda = density_active ? lambda_base * lambda_control : 0.0;
        const Real progress = iteration / static_cast<Real>(std::max(1, config.iterations));
        const Real fraction = density_active ? config.nonsmooth_step_fraction
                                             : config.nonsmooth_hpwl_step_fraction;
        Real learning_rate = fraction * ((db.xh - db.xl) + (db.yh - db.yl)) /
                             std::sqrt(1.0 + 4.0 * progress);
        if (refinement_active) learning_rate *= config.refinement_learning_rate_scale;
        const GlobalOptimizer active_optimizer = refinement_active
            ? config.refinement_optimizer : config.optimizer;
        BundleStats bundle_stats;
        if (config.enable_bundle && density_active && d.overflow <= config.bundle_start_overflow) {
            bundle_stats = apply_bundle(db, hpwl, learning_rate, iteration, config,
                                        cuts, last_cut, wgx, wgy);
        }

        Metrics current;
        current.iteration = iteration;
        current.exact_hpwl = hpwl;
        current.smooth_wirelength = hpwl;
        current.density_energy = d.energy;
        current.overflow = d.overflow;
        current.max_density = d.max_density;
        current.density_weight = lambda;
        current.step = learning_rate;
        if (d.overflow < best_overflow_metrics.overflow) {
            best_overflow_metrics = current;
            best_overflow = capture(db, fillers);
        }
        if (d.overflow <= config.stop_overflow &&
            (!result.have_feasible || hpwl < result.best_feasible_metrics.exact_hpwl)) {
            result.have_feasible = true;
            result.best_feasible_metrics = current;
            best_feasible = capture(db, fillers);
        }
        if (metrics) {
            metrics << iteration << ',' << hpwl << ',' << d.overflow << ','
                    << d.max_density << ',' << d.energy << ',' << lambda_base << ','
                    << lambda_control << ',' << lambda << ',' << learning_rate << ','
                    << bundle_stats.cuts << ',' << bundle_stats.newest_weight << ','
                    << bundle_stats.model_error << ',' << bundle_stats.cosine << ','
                    << (refinement_active ? 1 : 0) << ','
                    << global_optimizer_name(active_optimizer) << '\n';
        }
        if (iteration % config.log_every == 0 || iteration + 1 == config.iterations) {
            std::cout << "[Exact] iter=" << iteration << " hpwl=" << hpwl
                      << " overflow=" << d.overflow << " eta=" << lambda_control
                      << " lambda=" << lambda << " lr=" << learning_rate
                      << " cuts=" << bundle_stats.cuts
                      << " refine=" << (refinement_active ? 1 : 0) << '\n';
        }
        save_global_snapshot(db, config, iteration);

        std::vector<Real> gradient(2 * n, 0.0);
        std::size_t k = 0;
        for (int id : db.movable_ids) {
            const Node& node = db.nodes[id];
            const Real preconditioner = std::max<Real>(
                1.0, db.node_pin_weight[id] + lambda * node.area());
            gradient[k++] = (wgx[id] + lambda * dgx[id]) / preconditioner;
        }
        for (std::size_t i = 0; i < fillers.size(); ++i) {
            const Real preconditioner = std::max<Real>(
                1.0, lambda * fillers[i].width * fillers[i].height);
            gradient[k++] = lambda * fgx[i] / preconditioner;
        }
        k = n;
        for (int id : db.movable_ids) {
            const Node& node = db.nodes[id];
            const Real preconditioner = std::max<Real>(
                1.0, db.node_pin_weight[id] + lambda * node.area());
            gradient[k++] = (wgy[id] + lambda * dgy[id]) / preconditioner;
        }
        for (std::size_t i = 0; i < fillers.size(); ++i) {
            const Real preconditioner = std::max<Real>(
                1.0, lambda * fillers[i].width * fillers[i].height);
            gradient[k++] = lambda * fgy[i] / preconditioner;
        }

        std::vector<Real> positions = capture(db, fillers);
        const int age = refinement_active ? iteration - refinement_start + 1
                                          : iteration + 1;
        for (std::size_t i = 0; i < positions.size(); ++i) {
            positions[i] -= coordinate_step(active_optimizer, gradient[i], learning_rate,
                                             age, config, first[i], second[i],
                                             maximum_second[i]);
        }
        apply(db, fillers, positions);
        if (density_active) {
            ++density_step;
            if (refinement_active) {
                lambda_control = feasible_band_lambda(lambda_control, d.overflow, config);
            } else {
                lambda_control = update_lambda_control(
                    lambda_control, config.lambda_policy, density_step, hpwl, d.overflow,
                    previous_hpwl, initial_density_overflow, previous_overflow,
                    integral_error, config);
            }
        }
        previous_hpwl = hpwl;
    }

    if (result.have_feasible) {
        apply(db, fillers, best_feasible);
        std::cout << "[Exact] restored feasible iter="
                  << result.best_feasible_metrics.iteration << " hpwl="
                  << result.best_feasible_metrics.exact_hpwl << " overflow="
                  << result.best_feasible_metrics.overflow << '\n';
    } else {
        apply(db, fillers, best_overflow);
        std::cout << "[Exact] no feasible state; restored min-overflow iter="
                  << best_overflow_metrics.iteration << " overflow="
                  << best_overflow_metrics.overflow << '\n';
    }
    const DensityResult final_density = density.compute(
        db, fillers, nullptr, nullptr, nullptr, nullptr);
    result.final_metrics.exact_hpwl = exact_hpwl(db);
    result.final_metrics.smooth_wirelength = result.final_metrics.exact_hpwl;
    result.final_metrics.overflow = final_density.overflow;
    result.final_metrics.max_density = final_density.max_density;
    result.final_metrics.density_weight = lambda_base * lambda_control;
    result.wall_time_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return result;
}

}  // namespace dpcpp
