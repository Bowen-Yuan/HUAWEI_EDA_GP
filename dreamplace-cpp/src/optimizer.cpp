#include "optimizer.h"

#include "bookshelf.h"
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

struct DensityEvaluation {
    DensityResult metrics;
    std::vector<Real> node_gx;
    std::vector<Real> node_gy;
    std::vector<Real> filler_gx;
    std::vector<Real> filler_gy;
    bool have_gradient = false;
};

struct Evaluation {
    Metrics metrics;
    std::vector<Real> gradient;
    DensityEvaluation density;
};

void save_smooth_global_snapshot(const Database& db,
                                 const GlobalPlaceConfig& config,
                                 int iteration) {
    if (config.snapshot_every <= 0 || config.snapshot_dir.empty() ||
        iteration % config.snapshot_every != 0) return;
    std::ostringstream name;
    name << "global_" << std::setw(4) << std::setfill('0') << iteration
         << ".pl";
    write_bookshelf_pl(
        db, (std::filesystem::path(config.snapshot_dir) / name.str()).string());
}

std::vector<Real> capture(const Database& db, const std::vector<Filler>& fillers) {
    const std::size_t n = db.movable_ids.size() + fillers.size();
    std::vector<Real> result(2 * n);
    const std::size_t movable = db.movable_ids.size();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(movable); ++i) {
        const Node& node = db.nodes[db.movable_ids[i]];
        result[i] = node.x;
        result[n + i] = node.y;
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fillers.size()); ++i) {
        result[movable + i] = fillers[i].x;
        result[n + movable + i] = fillers[i].y;
    }
    return result;
}

void capture_into(const Database& db, const std::vector<Filler>& fillers,
                  std::vector<Real>& result) {
    const std::size_t n = db.movable_ids.size() + fillers.size();
    result.resize(2 * n);
    const std::size_t movable = db.movable_ids.size();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(movable); ++i) {
        const Node& node = db.nodes[db.movable_ids[i]];
        result[i] = node.x;
        result[n + i] = node.y;
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fillers.size()); ++i) {
        result[movable + i] = fillers[i].x;
        result[n + movable + i] = fillers[i].y;
    }
}

void apply(Database& db, std::vector<Filler>& fillers, const std::vector<Real>& values) {
    const std::size_t n = db.movable_ids.size() + fillers.size();
    if (values.size() != 2 * n) throw std::runtime_error("position vector size mismatch");
    const std::size_t movable = db.movable_ids.size();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(movable); ++i) {
        Node& node = db.nodes[db.movable_ids[i]];
        node.x = std::clamp(values[i], db.xl + 0.5 * node.width,
                            db.xh - 0.5 * node.width);
        node.y = std::clamp(values[n + i], db.yl + 0.5 * node.height,
                            db.yh - 0.5 * node.height);
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fillers.size()); ++i) {
        Filler& filler = fillers[i];
        filler.x = std::clamp(values[movable + i],
                              db.xl + 0.5 * filler.width,
                              db.xh - 0.5 * filler.width);
        filler.y = std::clamp(values[n + movable + i],
                              db.yl + 0.5 * filler.height,
                              db.yh - 0.5 * filler.height);
    }
}

Real l1_norm(const std::vector<Real>& x) {
    Real result = 0.0;
    for (Real value : x) result += std::abs(value);
    return result;
}

Real gamma_from_overflow(const ElectricDensity& density,
                         const GlobalPlaceConfig& config, Real overflow) {
    const Real base = config.gamma_scale *
        (density.bin_width() + density.bin_height());
    const Real exponent = (overflow - 0.1) * 20.0 / 9.0 - 1.0;
    return base * std::pow(10.0, std::clamp(exponent, -3.0, 3.0));
}

void evaluate_density(Database& db, std::vector<Filler>& fillers,
                      ElectricDensity& density, bool need_gradient,
                      DensityEvaluation& result) {
    result.have_gradient = need_gradient;
    result.metrics = density.compute(
        db, fillers,
        need_gradient ? &result.node_gx : nullptr,
        need_gradient ? &result.node_gy : nullptr,
        need_gradient ? &result.filler_gx : nullptr,
        need_gradient ? &result.filler_gy : nullptr);
}

void populate_evaluation(Database& db,
                         const std::vector<Filler>& fillers,
                         const GlobalPlaceConfig& config,
                         Real lambda, Real gamma, bool need_gradient,
                         Evaluation& result) {
    if (need_gradient && !result.density.have_gradient)
        throw std::runtime_error("cached density evaluation has no gradients");
    const DensityResult& d = result.density.metrics;
    static std::vector<Real> wire_gx;
    static std::vector<Real> wire_gy;
    const Real smooth = weighted_average_wirelength(
        db, gamma, config.degree_limit,
        need_gradient ? &wire_gx : nullptr,
        need_gradient ? &wire_gy : nullptr);
    result.metrics = Metrics{};
    result.metrics.exact_hpwl = exact_hpwl(db);
    result.metrics.smooth_wirelength = smooth;
    result.metrics.density_energy = d.energy;
    result.metrics.overflow = d.overflow;
    result.metrics.max_density = d.max_density;
    result.metrics.density_weight = lambda;
    result.metrics.gamma = gamma;
    if (!need_gradient) return;

    const std::size_t n = db.movable_ids.size() + fillers.size();
    result.gradient.resize(2 * n);
    const std::size_t movable = db.movable_ids.size();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(movable); ++i) {
        const int id = db.movable_ids[i];
        const Node& node = db.nodes[id];
        const Real preconditioner = std::max<Real>(
            1.0, db.node_pin_weight[id] + lambda * node.area());
        result.gradient[i] =
            (wire_gx[id] + lambda * result.density.node_gx[id]) / preconditioner;
        result.gradient[n + i] =
            (wire_gy[id] + lambda * result.density.node_gy[id]) / preconditioner;
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(fillers.size()); ++i) {
        const Real preconditioner = std::max<Real>(
            1.0, lambda * fillers[i].width * fillers[i].height);
        result.gradient[movable + i] =
            lambda * result.density.filler_gx[i] / preconditioner;
        result.gradient[n + movable + i] =
            lambda * result.density.filler_gy[i] / preconditioner;
    }
}

Evaluation evaluate(Database& db, std::vector<Filler>& fillers,
                    ElectricDensity& density, const GlobalPlaceConfig& config,
                    Real lambda, Real gamma, bool need_gradient) {
    Evaluation result;
    evaluate_density(db, fillers, density, need_gradient, result.density);
    populate_evaluation(db, fillers, config, lambda, gamma, need_gradient, result);
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
    case LambdaPolicy::BandDual: return "hybrid-band-dual";
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

    Real initial_ss = 0.0;
    Real initial_yy = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const Real s = v[i] - v_previous[i];
        const Real y = current.gradient[i] - previous_eval.gradient[i];
        initial_ss += s * s;
        initial_yy += y * y;
    }
    Real step = std::sqrt(initial_ss) /
                std::max<Real>(std::sqrt(initial_yy), 1.0e-30);
    Real acceleration = 1.0;
    Real previous_hpwl = current.metrics.exact_hpwl;
    std::vector<Real> best_feasible_position;
    std::vector<Real> best_overflow_position = v;
    Metrics best_overflow_metrics = current.metrics;
    std::vector<Real> u_next(v.size());
    std::vector<Real> v_next(v.size());

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
        save_smooth_global_snapshot(db, config, iteration);

        const Real next_acceleration =
            0.5 * (1.0 + std::sqrt(4.0 * acceleration * acceleration + 1.0));
        const Real momentum = (acceleration - 1.0) / next_acceleration;
        // DREAMPlace's BB implementation reevaluates both reference points
        // under the current lambda and gamma.  Reusing the previous
        // iteration's gradient here mixes two different objectives and makes
        // the curvature estimate meaningless as density weight changes.
        apply(db, fillers, v_previous);
        if (iteration > 0) {
            populate_evaluation(
                db, fillers, config, lambda, gamma, true, previous_eval);
            ++result.objective_evaluations;
        }
        apply(db, fillers, v);
        Real ss = 0.0;
        Real sy = 0.0;
        Real yy = 0.0;
        for (std::size_t i = 0; i < v.size(); ++i) {
            const Real s = v[i] - v_previous[i];
            const Real y = current.gradient[i] - previous_eval.gradient[i];
            ss += s * s;
            sy += s * y;
            yy += y * y;
        }
        const Real lip = std::sqrt(ss) /
                         std::max<Real>(std::sqrt(yy), 1.0e-30);
        const Real bb_short = sy / std::max<Real>(yy, 1.0e-30);
        if (std::isfinite(bb_short) && bb_short > 0.0) step = bb_short;
        else step = std::min(step, lip);
        // The upstream BB solver does not cap alpha by a fraction of the
        // placement region.  Keep only a very loose finite guard; boundary
        // projection handles an aggressive reference step.
        const Real max_step = 100.0 * std::max(db.xh - db.xl, db.yh - db.yl);
        step = std::clamp(step, 1.0e-12, max_step);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(v.size()); ++i) {
            u_next[i] = v[i] - step * current.gradient[i];
            v_next[i] = u_next[i] + momentum * (u_next[i] - u[i]);
        }
        v_previous.swap(v);
        std::swap(previous_eval, current);
        u.swap(u_next);
        apply(db, fillers, v_next);
        capture_into(db, fillers, v);
        acceleration = next_acceleration;

        lambda = update_density_weight(lambda,
            previous_eval.metrics.exact_hpwl - previous_hpwl,
            iteration, config.reference_hpwl_delta);
        previous_hpwl = previous_eval.metrics.exact_hpwl;
        const bool need_next_gradient = iteration + 1 < config.iterations;
        evaluate_density(
            db, fillers, density, need_next_gradient, current.density);
        gamma = gamma_from_overflow(
            density, config, current.density.metrics.overflow);
        populate_evaluation(
            db, fillers, config, lambda, gamma, need_next_gradient, current);
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
