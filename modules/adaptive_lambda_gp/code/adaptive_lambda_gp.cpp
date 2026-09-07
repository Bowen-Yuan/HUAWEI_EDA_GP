#include "adaptive_lambda_gp.hpp"

#include "epsilon_active/density.hpp"
#include "epsilon_active/hpwl.hpp"
#include "epsilon_active/optimizer.hpp"
#include "epsilon_active/step_policy.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace nsgp::adaptive {
namespace {

void clamp_movable(ea::Database& db) {
    for (const int id : db.movable_ids) {
        auto& node = db.nodes[id];
        node.x = std::clamp(node.x, db.xl + 0.5 * node.width, db.xh - 0.5 * node.width);
        node.y = std::clamp(node.y, db.yl + 0.5 * node.height, db.yh - 0.5 * node.height);
    }
}

Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
    Real value = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) value += a[i] * b[i];
    return value;
}

// RMS normalization: keeps lambda a dimensionless direction weight.
void rms_normalize(std::vector<Real>& x) {
    const Real norm = std::sqrt(std::max(1.0e-30, dot(x, x) /
                                               static_cast<Real>(std::max<std::size_t>(1, x.size()))));
    for (Real& value : x) value /= norm;
}

std::vector<Real> simplex_project(std::vector<Real> v) {
    std::vector<Real> u = v;
    std::sort(u.rbegin(), u.rend());
    Real total = 0.0, theta = 0.0;
    int rho = 0;
    for (std::size_t j = 0; j < u.size(); ++j) {
        total += u[j];
        const Real q = (total - 1.0) / static_cast<Real>(j + 1);
        if (u[j] > q) {
            rho = static_cast<int>(j) + 1;
            theta = q;
        }
    }
    if (!rho) return std::vector<Real>(v.size(), 1.0 / static_cast<Real>(v.size()));
    for (Real& x : v) x = std::max(static_cast<Real>(0.0), x - theta);
    return v;
}

// Minimum-norm convex combination on the simplex of the given directions.
std::vector<Real> min_norm_combination(const std::vector<std::vector<Real>>& dirs) {
    const std::size_t k = dirs.size();
    if (k == 0) throw std::runtime_error("empty direction ensemble");
    if (k == 1) return dirs[0];
    std::vector<Real> weights(k, 1.0 / static_cast<Real>(k)), gram(k * k);
    Real lipschitz = 1.0e-12;
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < k; ++j) {
            gram[i * k + j] = dot(dirs[i], dirs[j]);
            lipschitz = std::max(lipschitz, std::abs(gram[i * k + j]));
        }
    }
    for (int step = 0; step < 80; ++step) {
        std::vector<Real> grad(k, 0.0);
        for (std::size_t i = 0; i < k; ++i) {
            for (std::size_t j = 0; j < k; ++j) grad[i] += gram[i * k + j] * weights[j];
        }
        for (std::size_t i = 0; i < k; ++i) {
            weights[i] -= grad[i] / (lipschitz * static_cast<Real>(k));
        }
        weights = simplex_project(weights);
    }
    std::vector<Real> out(dirs[0].size(), 0.0);
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < out.size(); ++j) out[j] += weights[i] * dirs[i][j];
    }
    return out;
}

}  // namespace

AdaptiveLambdaRunStats run_adaptive_lambda_search(
    ea::Database& db, const int bins_x, const int bins_y, const Real target_density,
    const AdaptiveLambdaOptions& options, const bool telemetry) {
    options.validate();
    const std::size_t n = db.nodes.size();
    if (n == 0 || db.movable_ids.empty()) {
        throw std::runtime_error("adaptive_lambda_gp requires a non-empty layout");
    }

    ea::ExactHpwl hpwl(db);
    ea::ExactOverlapDensity density(db, bins_x, bins_y, target_density);
    const Real bin = std::min(density.bin_width(), density.bin_height());
    auto evaluate_exact = [&]() {
        const auto dm = density.evaluate(0.0, 1.0, nullptr, nullptr);
        FunnelMetrics metrics;
        metrics.hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        metrics.overflow = dm.overflow;
        return metrics;
    };

    auto optimizer = ea::make_optimizer(ea::parse_optimizer(options.optimizer),
                                        options.beta1, options.beta2,
                                        options.momentum, options.numerical_epsilon);
    optimizer->reset(2 * n);
    ea::StepControllerConfig step_config;
    step_config.learning_rate = options.learning_rate;
    step_config.maximum_delta_bins = options.maximum_delta_bins;
    step_config.learning_rate_min_ratio = options.learning_rate_min_ratio;
    step_config.trust_radius_bins = options.trust_radius_bins;
    step_config.radius_min_bins = options.radius_min_bins;
    step_config.radius_max_bins = options.radius_max_bins;
    step_config.grow_factor = options.grow_factor;
    step_config.shrink_factor = options.shrink_factor;
    step_config.grow_after_accepts = options.grow_after_accepts;
    auto step = ea::make_step_controller(ea::parse_step_policy(options.step_policy),
                                         step_config);
    step->reset(options.iterations, bin);

    AdaptiveLambdaRunStats stats;
    stats.initial = evaluate_exact();
    stats.best_any = stats.initial;
    // The input checkpoint is canonical feasible, so the best-feasible
    // fallback exists from iteration 0.
    stats.best_feasible = stats.initial;
    stats.maximum_observed_overflow = stats.initial.overflow;
    stats.lambda_initial = 0.0;

    std::vector<Real> best_feasible_positions(2 * n, 0.0);
    auto save_movable = [&](std::vector<Real>& out) {
        std::fill(out.begin(), out.end(), 0.0);
        for (const int id : db.movable_ids) {
            out[static_cast<std::size_t>(id)] = db.nodes[id].x;
            out[static_cast<std::size_t>(n) + static_cast<std::size_t>(id)] = db.nodes[id].y;
        }
    };
    auto restore_movable = [&](const std::vector<Real>& pos) {
        for (const int id : db.movable_ids) {
            db.nodes[id].x = pos[static_cast<std::size_t>(id)];
            db.nodes[id].y = pos[static_cast<std::size_t>(n) + static_cast<std::size_t>(id)];
        }
    };
    save_movable(best_feasible_positions);

    FunnelLambdaController lambda_controller(options.lambda);
    stats.lambda_initial = lambda_controller.lambda();
    stats.lambda_max_observed = lambda_controller.lambda();

    const Real final_ratio = options.lambda.final_overflow_ratio();
    const Real feasible_cap = final_ratio + AdaptiveLambdaOptions::kFeasibleToleranceRatio;
    FunnelMetrics current = stats.initial;
    int reject_streak = 0;
    bool final_lock_reset_done = false;
    bool been_above_final = current.overflow > feasible_cap;
    std::vector<Real> wire_x, wire_y, dense_x, dense_y, gradient, delta, previous;

    for (int it = 1; it <= options.iterations; ++it) {
        // Wire direction: exact active-face ensemble, RMS-normalized per
        // epsilon, combined by minimum norm on the simplex.
        std::vector<std::vector<Real>> wires;
        wires.reserve(options.epsilon_bin_scales.size());
        for (const Real scale : options.epsilon_bin_scales) {
            hpwl.evaluate(scale * bin, options.active_power, -1, &wire_x, &wire_y);
            std::vector<Real> direction(2 * n, 0.0);
            for (std::size_t i = 0; i < n; ++i) {
                direction[i] = wire_x[i];
                direction[n + i] = wire_y[i];
            }
            rms_normalize(direction);
            wires.push_back(std::move(direction));
            ++stats.objective_evaluations;
        }
        std::vector<Real> wire = min_norm_combination(wires);

        // Density direction: canonical exact overlap gradient, RMS
        // normalized independently before the lambda-weighted combination.
        density.evaluate(0.0, options.active_power, &dense_x, &dense_y);
        ++stats.objective_evaluations;
        std::vector<Real> dense(2 * n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            dense[i] = dense_x[i];
            dense[n + i] = dense_y[i];
        }
        rms_normalize(dense);

        gradient.resize(2 * n);
        const Real lambda = lambda_controller.lambda();
        for (std::size_t i = 0; i < 2 * n; ++i) {
            gradient[i] = wire[i] + lambda * dense[i];
        }
        rms_normalize(gradient);

        if (options.reset_on_final_lock && !final_lock_reset_done &&
            in_final_lock(options.lambda, it, options.iterations)) {
            optimizer->reset(2 * n);
            step->reset_streaks();
            ++stats.optimizer_resets;
            final_lock_reset_done = true;
        }

        const ea::StepDecision decision = step->propose(it);
        optimizer->compute_delta(gradient, decision.learning_rate,
                                 decision.maximum_delta, delta);

        previous.resize(2 * n);
        save_movable(previous);
        const Real corridor = corridor_ratio(options.lambda, it, options.iterations);
        FunnelMetrics candidate = current;
        bool accepted = false;
        int used_backtracks = 0;
        for (int backtrack = 0; backtrack < options.max_backtracks && !accepted; ++backtrack) {
            const Real alpha = std::ldexp(static_cast<Real>(1.0), -backtrack);
            for (const int id : db.movable_ids) {
                auto& node = db.nodes[id];
                node.x = previous[static_cast<std::size_t>(id)] -
                         alpha * delta[static_cast<std::size_t>(id)];
                node.y = previous[static_cast<std::size_t>(n) + static_cast<std::size_t>(id)] -
                         alpha * delta[static_cast<std::size_t>(n) + static_cast<std::size_t>(id)];
            }
            clamp_movable(db);
            candidate = evaluate_exact();
            stats.objective_evaluations += 2;
            used_backtracks = backtrack;
            accepted = funnel_accept(options.acceptance, corridor,
                                     stats.initial.hpwl, current, candidate) !=
                       FunnelDecision::Reject;
        }

        if (accepted) {
            current = candidate;
            ++stats.accepted;
            reject_streak = 0;
            if (candidate.hpwl < stats.best_any.hpwl) stats.best_any = candidate;
            if (candidate.overflow <= feasible_cap &&
                candidate.hpwl < stats.best_feasible.hpwl) {
                stats.best_feasible = candidate;
                save_movable(best_feasible_positions);
            }
            step->observe({true, used_backtracks});
        } else {
            restore_movable(previous);
            ++stats.rejected;
            ++reject_streak;
            step->observe({false, options.max_backtracks});
        }

        stats.maximum_observed_overflow =
            std::max(stats.maximum_observed_overflow, current.overflow);
        if (current.overflow > 0.10) ++stats.iterations_above_10_percent;
        if (been_above_final && current.overflow <= feasible_cap &&
            stats.first_return_to_feasible_iteration < 0) {
            stats.first_return_to_feasible_iteration = it;
        }
        been_above_final = been_above_final || current.overflow > feasible_cap;

        if (reject_streak >= options.reset_reject_streak) {
            optimizer->reset(2 * n);
            step->reset_streaks();
            reject_streak = 0;
            ++stats.optimizer_resets;
        }

        const Real lambda_before = lambda_controller.lambda();
        const Real lambda_ratio =
            lambda_controller.maybe_update(it, options.iterations, current.overflow);
        stats.lambda_max_observed =
            std::max(stats.lambda_max_observed, lambda_controller.lambda());
        if (lambda_ratio >= options.reset_lambda_ratio && lambda_ratio > 1.0) {
            optimizer->reset(2 * n);
            step->reset_streaks();
            ++stats.optimizer_resets;
        }
        (void)lambda_before;

        if (telemetry && options.verbose) {
            std::cout << "alambda iter=" << it << std::setprecision(14)
                      << " hpwl=" << current.hpwl
                      << " overflow_percent=" << current.overflow * 100.0
                      << " corridor_percent=" << corridor * 100.0
                      << " lambda=" << lambda_controller.lambda()
                      << " lr=" << decision.learning_rate
                      << " max_delta=" << decision.maximum_delta
                      << " accepted=" << (accepted ? 1 : 0)
                      << " backtracks=" << used_backtracks << '\n';
        }
    }

    // Stage-boundary contract: unconditionally fall back to the canonical
    // feasible best and re-audit it freshly.
    stats.last_before_restore = current;
    restore_movable(best_feasible_positions);
    clamp_movable(db);
    stats.final_selected = evaluate_exact();
    stats.lambda_final = lambda_controller.lambda();
    if (stats.final_selected.overflow > feasible_cap) {
        throw std::runtime_error(
            "adaptive_lambda_gp could not restore a canonical-feasible layout; "
            "the input checkpoint must be feasible under the final overflow target");
    }
    return stats;
}

}  // namespace nsgp::adaptive
