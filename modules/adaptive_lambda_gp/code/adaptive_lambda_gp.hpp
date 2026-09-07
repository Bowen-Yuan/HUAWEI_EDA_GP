#pragma once

// Pure search-policy logic for the adaptive_lambda_gp stage.  Everything in
// this header is dependency-free (types.hpp only) and header-only so the
// synthetic contract tests can exercise it without a benchmark.  The layout
// search itself lives in adaptive_lambda_gp.cpp.

#include "epsilon_active/types.hpp"

#include "epsilon_active/optimizer.hpp"
#include "epsilon_active/step_policy.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace nsgp::adaptive {

using Real = ea::Real;

// Overflow funnel + log-domain PI(D) lambda controller (plan V4 section 6).
struct FunnelLambdaConfig {
    bool fixed = false;  // ablation mode: keep lambda at "initial" forever
    Real initial = 0.25;
    Real minimum = 0.03;
    Real maximum = 16.0;
    int update_interval = 5;
    Real kp = 0.20;
    Real ki = 0.02;
    Real kd = 0.10;
    Real explore_overflow_percent = 15.0;
    Real final_overflow_percent = 7.0;
    Real hard_overflow_percent = 16.0;
    Real hold_fraction = 0.20;
    Real contract_end_fraction = 0.75;

    Real explore_overflow_ratio() const noexcept { return explore_overflow_percent * 0.01; }
    Real final_overflow_ratio() const noexcept { return final_overflow_percent * 0.01; }
    Real hard_overflow_ratio() const noexcept { return hard_overflow_percent * 0.01; }

    void validate() const {
        if (!(minimum > 0.0) || !(minimum <= initial) || !(initial <= maximum)) {
            throw std::invalid_argument("lambda bounds must satisfy 0 < min <= initial <= max");
        }
        if (update_interval < 1) throw std::invalid_argument("lambda update_interval must be >= 1");
        if (!(kp >= 0.0) || !(ki >= 0.0) || !(kd >= 0.0)) {
            throw std::invalid_argument("lambda gains must be non-negative");
        }
        const Real explore = explore_overflow_ratio();
        const Real final = final_overflow_ratio();
        const Real hard = hard_overflow_ratio();
        if (!(final > 0.0) || !(explore >= final) || !(hard >= explore)) {
            throw std::invalid_argument("overflow percentages must satisfy 0 < final <= explore <= hard");
        }
        if (!(hold_fraction >= 0.0) || !(hold_fraction < 1.0) ||
            !(contract_end_fraction > hold_fraction) || !(contract_end_fraction <= 1.0)) {
            throw std::invalid_argument("funnel fractions must satisfy 0 <= hold < end <= 1");
        }
    }
};

// Corridor height (overflow ratio) at a given iteration.  Phase A holds the
// explore height, phase B contracts with a smoothstep, phase C locks at the
// final target.  This is an acceptance funnel only; it never redefines the
// canonical density evaluator.
inline Real corridor_ratio(const FunnelLambdaConfig& config, int iteration,
                           int total_iterations) {
    const Real p = total_iterations > 0
        ? static_cast<Real>(iteration) / static_cast<Real>(total_iterations)
        : static_cast<Real>(1.0);
    const Real explore = config.explore_overflow_ratio();
    const Real final = config.final_overflow_ratio();
    if (p <= config.hold_fraction) return explore;
    if (p >= config.contract_end_fraction) return final;
    const Real u = (p - config.hold_fraction) /
                   (config.contract_end_fraction - config.hold_fraction);
    const Real s = u * u * (3.0 - 2.0 * u);
    return final + (explore - final) * (1.0 - s);
}

inline bool in_final_lock(const FunnelLambdaConfig& config, int iteration,
                          int total_iterations) {
    const Real p = total_iterations > 0
        ? static_cast<Real>(iteration) / static_cast<Real>(total_iterations)
        : static_cast<Real>(1.0);
    return p >= config.contract_end_fraction;
}

class FunnelLambdaController {
public:
    explicit FunnelLambdaController(FunnelLambdaConfig config)
        : config_(config), lambda_(config.initial) {
        config_.validate();
    }

    void reset() {
        lambda_ = config_.initial;
        integral_ = 0.0;
        has_previous_ = false;
        previous_overflow_ = 0.0;
    }

    Real lambda() const noexcept { return lambda_; }
    Real integral() const noexcept { return integral_; }
    Real last_error() const noexcept { return last_error_; }
    Real last_trend() const noexcept { return last_trend_; }

    // Applies at most one lambda update (every update_interval iterations).
    // Returns lambda_new / lambda_old for the caller's reset rule; returns
    // exactly 1.0 on iterations without an update.
    Real maybe_update(int iteration, int total_iterations,
                      Real current_overflow_ratio) {
        if (config_.fixed) return 1.0;
        if (iteration % config_.update_interval != 0) return 1.0;
        const Real corridor = corridor_ratio(config_, iteration, total_iterations);
        constexpr Real kPercent = 0.01;  // 1% expressed as a ratio
        const Real error = (current_overflow_ratio - corridor) / kPercent;
        const Real trend = has_previous_
            ? (current_overflow_ratio - previous_overflow_) / kPercent
            : 0.0;
        integral_ = std::clamp(0.8 * integral_ + error, -4.0, 4.0);
        const Real log_step = std::clamp(
            config_.kp * error + config_.ki * integral_ + config_.kd * trend,
            -0.35, 0.55);
        const Real before = lambda_;
        lambda_ *= std::exp(log_step);
        if (in_final_lock(config_, iteration, total_iterations) &&
            current_overflow_ratio > config_.final_overflow_ratio()) {
            const Real extra = std::clamp(
                0.10 * (current_overflow_ratio - config_.final_overflow_ratio()) / kPercent,
                0.0, 0.25);
            lambda_ *= std::exp(extra);
        }
        lambda_ = std::clamp(lambda_, config_.minimum, config_.maximum);
        previous_overflow_ = current_overflow_ratio;
        has_previous_ = true;
        last_error_ = error;
        last_trend_ = trend;
        return lambda_ / before;
    }

private:
    FunnelLambdaConfig config_;
    Real lambda_;
    Real integral_ = 0.0;
    Real previous_overflow_ = 0.0;
    bool has_previous_ = false;
    Real last_error_ = 0.0;
    Real last_trend_ = 0.0;
};

// Exact funnel acceptance (plan V4 section 7).  Operates on canonical exact
// metrics only; the caller performs the fresh audits.
struct FunnelAcceptanceConfig {
    Real final_overflow_ratio = 0.07;
    Real hard_overflow_ratio = 0.16;
    Real hpwl_tolerance_relative = 1.0e-12;
    Real overflow_improvement_epsilon_ratio = 2.0e-5;  // 0.002%
    Real recovery_step_hpwl_budget = 5.0e-4;           // 0.05%
    Real recovery_total_hpwl_budget = 3.0e-3;          // 0.30%

    void validate() const {
        if (!(final_overflow_ratio > 0.0) || !(hard_overflow_ratio > final_overflow_ratio)) {
            throw std::invalid_argument("acceptance needs 0 < final < hard overflow");
        }
        if (!(hpwl_tolerance_relative >= 0.0)) {
            throw std::invalid_argument("hpwl tolerance must be non-negative");
        }
        if (!(overflow_improvement_epsilon_ratio > 0.0)) {
            throw std::invalid_argument("overflow improvement epsilon must be positive");
        }
        if (!(recovery_step_hpwl_budget >= 0.0) || !(recovery_total_hpwl_budget >= 0.0)) {
            throw std::invalid_argument("recovery hpwl budgets must be non-negative");
        }
    }
};

struct FunnelMetrics {
    Real hpwl = 0.0;
    Real overflow = 0.0;
};

enum class FunnelDecision {
    Reject,
    ExploreAccept,
    RecoveryAccept
};

inline FunnelDecision funnel_accept(const FunnelAcceptanceConfig& config,
                                    Real corridor, Real initial_hpwl,
                                    const FunnelMetrics& current,
                                    const FunnelMetrics& candidate) {
    if (!std::isfinite(candidate.hpwl) || !std::isfinite(candidate.overflow)) {
        return FunnelDecision::Reject;
    }
    if (candidate.overflow > config.hard_overflow_ratio) {
        return FunnelDecision::Reject;
    }
    const Real tolerance =
        config.hpwl_tolerance_relative * std::max<Real>(1.0, std::abs(current.hpwl));
    if (current.overflow <= corridor) {
        // Exploration: stay inside the funnel and strictly reduce HPWL.
        if (candidate.overflow <= corridor &&
            candidate.hpwl < current.hpwl - tolerance) {
            return FunnelDecision::ExploreAccept;
        }
        return FunnelDecision::Reject;
    }
    // Recovery: strictly reduce overflow under bounded HPWL cost.
    if (candidate.overflow >= current.overflow - config.overflow_improvement_epsilon_ratio) {
        return FunnelDecision::Reject;
    }
    if (candidate.hpwl > current.hpwl * (1.0 + config.recovery_step_hpwl_budget)) {
        return FunnelDecision::Reject;
    }
    if (candidate.hpwl > initial_hpwl * (1.0 + config.recovery_total_hpwl_budget)) {
        return FunnelDecision::Reject;
    }
    return FunnelDecision::RecoveryAccept;
}

struct AdaptiveLambdaRunStats {
    int accepted = 0;
    int rejected = 0;
    int objective_evaluations = 0;
    int optimizer_resets = 0;
    int iterations_above_10_percent = 0;
    int first_return_to_feasible_iteration = -1;
    Real maximum_observed_overflow = 0.0;
    Real lambda_initial = 0.0;
    Real lambda_max_observed = 0.0;
    Real lambda_final = 0.0;
    FunnelMetrics initial;
    FunnelMetrics best_any;
    FunnelMetrics best_feasible;
    FunnelMetrics last_before_restore;
    FunnelMetrics final_selected;
};

struct AdaptiveLambdaOptions {
    int iterations = 200;
    bool verbose = false;

    // direction
    bool active_ensemble = true;
    std::vector<Real> epsilon_bin_scales = {0.0, 0.25, 0.5, 1.0};
    Real active_power = 1.0;

    // optimizer
    std::string optimizer = "adam";
    Real beta1 = 0.9;
    Real beta2 = 0.999;
    Real momentum = 0.9;
    Real numerical_epsilon = 1.0e-8;
    Real learning_rate = 0.02;

    // step policy
    std::string step_policy = "trust";
    Real maximum_delta_bins = 0.5;
    Real learning_rate_min_ratio = 0.10;
    Real trust_radius_bins = 1.0;
    Real radius_min_bins = 0.02;
    Real radius_max_bins = 2.0;
    Real grow_factor = 1.25;
    Real shrink_factor = 0.5;
    int grow_after_accepts = 3;
    int max_backtracks = 9;

    // lambda funnel
    FunnelLambdaConfig lambda;

    // acceptance
    FunnelAcceptanceConfig acceptance;

    // optimizer reset rules
    int reset_reject_streak = 4;
    Real reset_lambda_ratio = 4.0;
    bool reset_on_final_lock = true;

    // Numerical tolerance for the stage-boundary feasibility contract.
    static constexpr Real kFeasibleToleranceRatio = 1.0e-9;

    void validate() const {
        if (iterations < 1) throw std::invalid_argument("iterations must be >= 1");
        ea::parse_optimizer(optimizer);
        (void)ea::parse_step_policy(step_policy);
        if (!(learning_rate > 0.0) || !(numerical_epsilon > 0.0)) {
            throw std::invalid_argument("learning rate and numerical epsilon must be positive");
        }
        if (!(beta1 >= 0.0) || !(beta1 < 1.0) || !(beta2 >= 0.0) || !(beta2 < 1.0) ||
            !(momentum >= 0.0) || !(momentum < 1.0)) {
            throw std::invalid_argument("invalid optimizer hyperparameters");
        }
        if (!(maximum_delta_bins > 0.0) || !(trust_radius_bins > 0.0)) {
            throw std::invalid_argument("step sizes must be positive");
        }
        if (!(learning_rate_min_ratio > 0.0) || !(learning_rate_min_ratio < 1.0)) {
            throw std::invalid_argument("learning_rate_min_ratio must be in (0, 1)");
        }
        if (!(radius_min_bins > 0.0) || !(radius_max_bins >= radius_min_bins) ||
            !(grow_factor > 1.0) || !(shrink_factor > 0.0) || !(shrink_factor < 1.0) ||
            grow_after_accepts < 1) {
            throw std::invalid_argument("invalid trust radius configuration");
        }
        if (max_backtracks < 0 || max_backtracks > 60) {
            throw std::invalid_argument("max_backtracks must be in [0, 60]");
        }
        if (epsilon_bin_scales.empty()) {
            throw std::invalid_argument("epsilon_bin_scales must not be empty");
        }
        for (const Real scale : epsilon_bin_scales) {
            if (!std::isfinite(scale) || !(scale >= 0.0)) {
                throw std::invalid_argument("epsilon_bin_scales entries must be >= 0");
            }
        }
        if (!(active_power > 0.0)) throw std::invalid_argument("active_power must be positive");
        if (reset_reject_streak < 1) throw std::invalid_argument("reset_reject_streak must be >= 1");
        if (!(reset_lambda_ratio >= 1.0)) throw std::invalid_argument("reset_lambda_ratio must be >= 1");
        lambda.validate();
        acceptance.validate();
    }
};

// Runs the in-memory search on db.  On return the movable layout is the
// canonical-feasible best (restored + clamped + fresh-audited).  telemetry,
// when non-null and options.verbose, receives one key=value line per
// iteration on stdout; no fourth experiment file is produced.
AdaptiveLambdaRunStats run_adaptive_lambda_search(
    ea::Database& db, int bins_x, int bins_y, Real target_density,
    const AdaptiveLambdaOptions& options, bool telemetry);

}  // namespace nsgp::adaptive
