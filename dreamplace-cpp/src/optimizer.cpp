#include "optimizer.h"

#include "nonsmooth.h"
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

namespace dpcpp {
namespace {

struct Evaluation {
    Metrics metrics;
    std::vector<Real> gradient;
};

std::vector<Real> capture(const Database& db, const std::vector<Filler>& fillers) {
    const std::size_t n = db.movable_ids.size() + fillers.size();
    std::vector<Real> result(2 * n);
    std::size_t k = 0;
    for (int id : db.movable_ids) result[k++] = db.nodes[id].x;
    for (const Filler& filler : fillers) result[k++] = filler.x;
    k = n;
    for (int id : db.movable_ids) result[k++] = db.nodes[id].y;
    for (const Filler& filler : fillers) result[k++] = filler.y;
    return result;
}

void apply(Database& db, std::vector<Filler>& fillers, const std::vector<Real>& values) {
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

Real l1_norm(const std::vector<Real>& x) {
    Real result = 0.0;
    for (Real value : x) result += std::abs(value);
    return result;
}

Real l2_norm(const std::vector<Real>& x) {
    Real result = 0.0;
    for (Real value : x) result += value * value;
    return std::sqrt(result);
}

Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) result += a[i] * b[i];
    return result;
}

std::vector<Real> difference(const std::vector<Real>& a,
                             const std::vector<Real>& b) {
    std::vector<Real> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) result[i] = a[i] - b[i];
    return result;
}

Real gamma_from_overflow(const ElectricDensity& density,
                         const GlobalPlaceConfig& config, Real overflow) {
    const Real base = config.gamma_scale *
        (density.bin_width() + density.bin_height());
    const Real exponent = (overflow - 0.1) * 20.0 / 9.0 - 1.0;
    return base * std::pow(10.0, std::clamp(exponent, -3.0, 3.0));
}

Evaluation evaluate(Database& db, std::vector<Filler>& fillers,
                    ElectricDensity& density, const GlobalPlaceConfig& config,
                    Real lambda, Real gamma, bool need_gradient) {
    std::vector<Real> density_gx, density_gy, filler_gx, filler_gy;
    const DensityResult d = density.compute(
        db, fillers,
        need_gradient ? &density_gx : nullptr,
        need_gradient ? &density_gy : nullptr,
        need_gradient ? &filler_gx : nullptr,
        need_gradient ? &filler_gy : nullptr);
    std::vector<Real> wire_gx, wire_gy;
    const Real smooth = weighted_average_wirelength(
        db, gamma, config.degree_limit,
        need_gradient ? &wire_gx : nullptr,
        need_gradient ? &wire_gy : nullptr);
    Evaluation result;
    result.metrics.exact_hpwl = exact_hpwl(db);
    result.metrics.smooth_wirelength = smooth;
    result.metrics.density_energy = d.energy;
    result.metrics.overflow = d.overflow;
    result.metrics.max_density = d.max_density;
    result.metrics.density_weight = lambda;
    result.metrics.gamma = gamma;
    if (!need_gradient) return result;

    const std::size_t n = db.movable_ids.size() + fillers.size();
    result.gradient.assign(2 * n, 0.0);
    std::size_t k = 0;
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real preconditioner = std::max<Real>(
            1.0, db.node_pin_weight[id] + lambda * node.area());
        result.gradient[k++] = (wire_gx[id] + lambda * density_gx[id]) /
                               preconditioner;
    }
    for (std::size_t i = 0; i < fillers.size(); ++i) {
        const Real preconditioner = std::max<Real>(
            1.0, lambda * fillers[i].width * fillers[i].height);
        result.gradient[k++] = lambda * filler_gx[i] / preconditioner;
    }
    k = n;
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real preconditioner = std::max<Real>(
            1.0, db.node_pin_weight[id] + lambda * node.area());
        result.gradient[k++] = (wire_gy[id] + lambda * density_gy[id]) /
                               preconditioner;
    }
    for (std::size_t i = 0; i < fillers.size(); ++i) {
        const Real preconditioner = std::max<Real>(
            1.0, lambda * fillers[i].width * fillers[i].height);
        result.gradient[k++] = lambda * filler_gy[i] / preconditioner;
    }
    return result;
}

Real initialize_density_weight(Database& db, std::vector<Filler>& fillers,
                               ElectricDensity& density,
                               const GlobalPlaceConfig& config, Real gamma) {
    std::vector<Real> density_gx, density_gy, filler_gx, filler_gy;
    density.compute(db, fillers, &density_gx, &density_gy, &filler_gx, &filler_gy);
    std::vector<Real> wire_gx, wire_gy;
    weighted_average_wirelength(db, gamma, config.degree_limit, &wire_gx, &wire_gy);
    Real wire_l1 = 0.0;
    Real density_l1 = 0.0;
    for (int id : db.movable_ids) {
        wire_l1 += std::abs(wire_gx[id]) + std::abs(wire_gy[id]);
        density_l1 += std::abs(density_gx[id]) + std::abs(density_gy[id]);
    }
    density_l1 += l1_norm(filler_gx) + l1_norm(filler_gy);
    const Real lambda = config.density_weight_scale * wire_l1 /
                        std::max<Real>(density_l1, 1.0e-30);
    std::cout << "[Lambda] wire_grad_l1=" << wire_l1
              << " density_grad_l1=" << density_l1
              << " initial=" << lambda << '\n';
    return lambda;
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

Real update_density_weight(Real lambda, Real delta_hpwl, int iteration,
                           Real reference_hpwl_delta) {
    constexpr Real lower = 0.95;
    constexpr Real upper = 1.05;
    Real mu = 1.0;
    if (delta_hpwl < 0.0) {
        mu = upper * std::max(std::pow(0.9999, iteration), 0.98);
    } else {
        mu = std::max(lower,
            std::pow(upper, 1.0 - delta_hpwl / reference_hpwl_delta));
    }
    return lambda * mu;
}

}  // namespace

const char* wirelength_model_name(WirelengthModel model) {
    switch (model) {
    case WirelengthModel::WeightedAverage: return "weighted-average";
    case WirelengthModel::ExactHpwl: return "exact-hpwl";
    }
    return "unknown";
}

const char* global_optimizer_name(GlobalOptimizer optimizer) {
    switch (optimizer) {
    case GlobalOptimizer::DreamplaceNesterov: return "dreamplace";
    case GlobalOptimizer::HeavyBall: return "heavy-ball";
    case GlobalOptimizer::Adam: return "adam";
    case GlobalOptimizer::AMSGrad: return "amsgrad";
    case GlobalOptimizer::AdaGrad: return "adagrad";
    }
    return "unknown";
}

const char* lambda_policy_name(LambdaPolicy policy) {
    switch (policy) {
    case LambdaPolicy::Dreamplace: return "dreamplace";
    case LambdaPolicy::Trajectory: return "trajectory";
    case LambdaPolicy::Ratio: return "ratio";
    }
    return "unknown";
}

GlobalPlaceResult global_place(Database& db, std::vector<Filler>& fillers,
                               const GlobalPlaceConfig& config) {
    if (config.wirelength_model == WirelengthModel::ExactHpwl) {
        if (config.optimizer == GlobalOptimizer::DreamplaceNesterov) {
            throw std::runtime_error(
                "exact HPWL is incompatible with BB-Nesterov; select "
                "heavy-ball, adam, amsgrad, or adagrad");
        }
        return nonsmooth_global_place(db, fillers, config);
    }
    if (config.optimizer != GlobalOptimizer::DreamplaceNesterov) {
        throw std::runtime_error(
            "the smooth weighted-average model currently requires the "
            "DREAMPlace BB-Nesterov optimizer");
    }
    const auto start_time = std::chrono::steady_clock::now();
    ElectricDensity density(db, config.bins_x, config.bins_y,
                            config.target_density);
    add_gp_noise(db, fillers, config.seed, config.gp_noise_ratio);

    DensityResult initial_density = density.compute(
        db, fillers, nullptr, nullptr, nullptr, nullptr);
    Real gamma = gamma_from_overflow(density, config, initial_density.overflow);
    Real lambda = initialize_density_weight(db, fillers, density, config, gamma);

    std::ofstream metrics_file;
    if (!config.metrics_path.empty()) {
        std::filesystem::create_directories(
            std::filesystem::path(config.metrics_path).parent_path());
        metrics_file.open(config.metrics_path);
        metrics_file << "iteration,exact_hpwl,wa_wirelength,overflow,max_density,"
                        "density_energy,lambda,gamma,step\n";
        metrics_file << std::setprecision(12);
    }

    GlobalPlaceResult result;
    std::vector<Real> v = capture(db, fillers);
    Evaluation current = evaluate(db, fillers, density, config, lambda, gamma, true);
    ++result.objective_evaluations;
    std::vector<Real> u = v;
    std::vector<Real> v_previous = v;
    for (std::size_t i = 0; i < v_previous.size(); ++i) {
        v_previous[i] -= config.initial_learning_rate * current.gradient[i];
    }
    apply(db, fillers, v_previous);
    v_previous = capture(db, fillers);
    Evaluation previous_eval = evaluate(db, fillers, density, config, lambda, gamma, true);
    ++result.objective_evaluations;
    apply(db, fillers, v);

    const auto initial_s = difference(v, v_previous);
    const auto initial_y = difference(current.gradient, previous_eval.gradient);
    Real step = l2_norm(initial_s) /
                std::max<Real>(l2_norm(initial_y), 1.0e-30);
    Real acceleration = 1.0;
    Real previous_hpwl = current.metrics.exact_hpwl;
    std::vector<Real> best_feasible_position;
    std::vector<Real> best_overflow_position = v;
    Metrics best_overflow_metrics = current.metrics;

    for (int iteration = 0; iteration < config.iterations; ++iteration) {
        current.metrics.iteration = iteration;
        current.metrics.step = step;
        if (current.metrics.overflow < best_overflow_metrics.overflow) {
            best_overflow_metrics = current.metrics;
            best_overflow_position = v;
        }
        if (current.metrics.overflow <= config.stop_overflow &&
            (!result.have_feasible ||
             current.metrics.exact_hpwl < result.best_feasible_metrics.exact_hpwl)) {
            result.have_feasible = true;
            result.best_feasible_metrics = current.metrics;
            best_feasible_position = v;
        }
        if (metrics_file) {
            metrics_file << iteration << ',' << current.metrics.exact_hpwl << ','
                         << current.metrics.smooth_wirelength << ','
                         << current.metrics.overflow << ','
                         << current.metrics.max_density << ','
                         << current.metrics.density_energy << ',' << lambda << ','
                         << gamma << ',' << step << '\n';
        }
        if (iteration % config.log_every == 0 || iteration + 1 == config.iterations) {
            std::cout << "[GP] iter=" << iteration
                      << " hpwl=" << current.metrics.exact_hpwl
                      << " overflow=" << current.metrics.overflow
                      << " max_density=" << current.metrics.max_density
                      << " lambda=" << lambda
                      << " gamma=" << gamma
                      << " step=" << step << '\n';
        }

        const Real next_acceleration =
            0.5 * (1.0 + std::sqrt(4.0 * acceleration * acceleration + 1.0));
        const Real momentum = (acceleration - 1.0) / next_acceleration;
        // DREAMPlace's BB implementation reevaluates both reference points
        // under the current lambda and gamma.  Reusing the previous
        // iteration's gradient here mixes two different objectives and makes
        // the curvature estimate meaningless as density weight changes.
        apply(db, fillers, v_previous);
        previous_eval = evaluate(db, fillers, density, config, lambda, gamma, true);
        ++result.objective_evaluations;
        apply(db, fillers, v);
        current = evaluate(db, fillers, density, config, lambda, gamma, true);
        ++result.objective_evaluations;
        const std::vector<Real> s = difference(v, v_previous);
        const std::vector<Real> y = difference(current.gradient,
                                               previous_eval.gradient);
        const Real sy = dot(s, y);
        const Real yy = dot(y, y);
        const Real lip = l2_norm(s) / std::max<Real>(l2_norm(y), 1.0e-30);
        const Real bb_short = sy / std::max<Real>(yy, 1.0e-30);
        if (std::isfinite(bb_short) && bb_short > 0.0) step = bb_short;
        else step = std::min(step, lip);
        // The upstream BB solver does not cap alpha by a fraction of the
        // placement region.  Keep only a very loose finite guard; boundary
        // projection handles an aggressive reference step.
        const Real max_step = 100.0 * std::max(db.xh - db.xl, db.yh - db.yl);
        step = std::clamp(step, 1.0e-12, max_step);

        std::vector<Real> u_next(v.size());
        std::vector<Real> v_next(v.size());
        for (std::size_t i = 0; i < v.size(); ++i) {
            u_next[i] = v[i] - step * current.gradient[i];
            v_next[i] = u_next[i] + momentum * (u_next[i] - u[i]);
        }
        v_previous = v;
        previous_eval = std::move(current);
        u = std::move(u_next);
        apply(db, fillers, v_next);
        v = capture(db, fillers);
        acceleration = next_acceleration;

        lambda = update_density_weight(lambda,
            previous_eval.metrics.exact_hpwl - previous_hpwl,
            iteration, config.reference_hpwl_delta);
        previous_hpwl = previous_eval.metrics.exact_hpwl;
        DensityResult d = density.compute(db, fillers, nullptr, nullptr, nullptr, nullptr);
        gamma = gamma_from_overflow(density, config, d.overflow);
        // The next loop only logs/checkpoints these metrics before step_bb
        // reevaluates the reference point under the current objective.  Do
        // not build a third, immediately-discarded gradient here.
        current = evaluate(db, fillers, density, config, lambda, gamma, false);
        ++result.objective_evaluations;
    }

    if (result.have_feasible) {
        apply(db, fillers, best_feasible_position);
        std::cout << "[GP] restored best feasible iter="
                  << result.best_feasible_metrics.iteration
                  << " hpwl=" << result.best_feasible_metrics.exact_hpwl
                  << " overflow=" << result.best_feasible_metrics.overflow << '\n';
    } else {
        apply(db, fillers, best_overflow_position);
        std::cout << "[GP] no <= " << config.stop_overflow
                  << " state; restored minimum-overflow iter="
                  << best_overflow_metrics.iteration
                  << " overflow=" << best_overflow_metrics.overflow << '\n';
    }
    DensityResult final_density = density.compute(
        db, fillers, nullptr, nullptr, nullptr, nullptr);
    result.final_metrics.exact_hpwl = exact_hpwl(db);
    result.final_metrics.overflow = final_density.overflow;
    result.final_metrics.max_density = final_density.max_density;
    result.final_metrics.density_weight = lambda;
    result.final_metrics.gamma = gamma;
    result.wall_time_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_time).count();
    return result;
}

}  // namespace dpcpp
