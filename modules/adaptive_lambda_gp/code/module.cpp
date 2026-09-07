#include "adaptive_lambda_gp.hpp"
#include "microkernel.hpp"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace {

using nsgp::adaptive::AdaptiveLambdaOptions;
using ea::Real;
using Json = nlohmann::json;

Real percent_to_ratio(const Json& node, const char* key, Real fallback) {
    return node.value(key, fallback) * 0.01;
}

AdaptiveLambdaOptions parse_options(const Json& config) {
    AdaptiveLambdaOptions options;
    options.iterations = config.value("iterations", options.iterations);
    options.verbose = config.value("verbose", options.verbose);

    if (config.contains("direction")) {
        const Json& direction = config.at("direction");
        options.active_ensemble = direction.value("active_ensemble", options.active_ensemble);
        if (direction.contains("epsilon_bin_scales")) {
            options.epsilon_bin_scales =
                direction.at("epsilon_bin_scales").get<std::vector<Real>>();
        }
        options.active_power = direction.value("active_power", options.active_power);
    }

    if (config.contains("optimizer")) {
        const Json& optimizer = config.at("optimizer");
        options.optimizer = optimizer.value("name", options.optimizer);
        options.learning_rate = optimizer.value("learning_rate", options.learning_rate);
        options.beta1 = optimizer.value("beta1", options.beta1);
        options.beta2 = optimizer.value("beta2", options.beta2);
        options.momentum = optimizer.value("momentum", options.momentum);
        options.numerical_epsilon = optimizer.value("numerical_epsilon", options.numerical_epsilon);
    }

    if (config.contains("step_policy")) {
        const Json& step = config.at("step_policy");
        options.step_policy = step.value("name", options.step_policy);
        options.maximum_delta_bins = step.value("maximum_delta_bins", options.maximum_delta_bins);
        options.learning_rate_min_ratio =
            step.value("learning_rate_min_ratio", options.learning_rate_min_ratio);
        options.trust_radius_bins = step.value("trust_radius_bins", options.trust_radius_bins);
        options.radius_min_bins = step.value("radius_min_bins", options.radius_min_bins);
        options.radius_max_bins = step.value("radius_max_bins", options.radius_max_bins);
        options.grow_factor = step.value("grow_factor", options.grow_factor);
        options.shrink_factor = step.value("shrink_factor", options.shrink_factor);
        options.grow_after_accepts = step.value("grow_after_accepts", options.grow_after_accepts);
        options.max_backtracks = step.value("max_backtracks", options.max_backtracks);
    }

    if (config.contains("lambda_policy")) {
        const Json& policy = config.at("lambda_policy");
        const std::string name = policy.value("name", std::string("funnel_pid"));
        if (name != "funnel_pid" && name != "fixed") {
            throw std::runtime_error("lambda_policy.name must be funnel_pid or fixed");
        }
        options.lambda.fixed = (name == "fixed");
        options.lambda.initial = policy.value("initial", options.lambda.initial);
        options.lambda.minimum = policy.value("minimum", options.lambda.minimum);
        options.lambda.maximum = policy.value("maximum", options.lambda.maximum);
        options.lambda.update_interval =
            policy.value("update_interval", options.lambda.update_interval);
        options.lambda.kp = policy.value("kp", options.lambda.kp);
        options.lambda.ki = policy.value("ki", options.lambda.ki);
        options.lambda.kd = policy.value("kd", options.lambda.kd);
        options.lambda.explore_overflow_percent = policy.value(
            "explore_overflow_percent", options.lambda.explore_overflow_percent);
        options.lambda.final_overflow_percent = policy.value(
            "final_overflow_percent", options.lambda.final_overflow_percent);
        options.lambda.hard_overflow_percent = policy.value(
            "hard_overflow_percent", options.lambda.hard_overflow_percent);
        options.lambda.hold_fraction = policy.value("hold_fraction", options.lambda.hold_fraction);
        options.lambda.contract_end_fraction = policy.value(
            "contract_end_fraction", options.lambda.contract_end_fraction);
    }

    if (config.contains("acceptance")) {
        const Json& acceptance = config.at("acceptance");
        const std::string policy = acceptance.value("policy", std::string("exact_funnel"));
        if (policy != "exact_funnel") {
            throw std::runtime_error("acceptance.policy must be exact_funnel");
        }
        options.acceptance.hpwl_tolerance_relative = acceptance.value(
            "hpwl_tolerance_relative", options.acceptance.hpwl_tolerance_relative);
        options.acceptance.overflow_improvement_epsilon_ratio = percent_to_ratio(
            acceptance, "overflow_improvement_epsilon_percent",
            options.acceptance.overflow_improvement_epsilon_ratio * 100.0);
        options.acceptance.recovery_step_hpwl_budget = percent_to_ratio(
            acceptance, "recovery_step_hpwl_budget_percent",
            options.acceptance.recovery_step_hpwl_budget * 100.0);
        options.acceptance.recovery_total_hpwl_budget = percent_to_ratio(
            acceptance, "recovery_total_hpwl_budget_percent",
            options.acceptance.recovery_total_hpwl_budget * 100.0);
    }

    if (config.contains("reset")) {
        const Json& reset = config.at("reset");
        options.reset_reject_streak = reset.value("reject_streak", options.reset_reject_streak);
        options.reset_lambda_ratio = reset.value("lambda_ratio", options.reset_lambda_ratio);
        options.reset_on_final_lock =
            reset.value("reset_on_final_lock", options.reset_on_final_lock);
    }

    options.validate();
    return options;
}

void print_metrics(const char* label, const nsgp::adaptive::FunnelMetrics& metrics) {
    std::cout << "[adaptive_lambda_gp] " << label << std::setprecision(14)
              << " hpwl=" << metrics.hpwl
              << " overflow_percent=" << metrics.overflow * 100.0 << '\n';
}

}  // namespace

namespace nsgp::modules {

void register_adaptive_lambda_gp(ModuleRegistry& registry) {
    registry.add("adaptive_lambda_gp", [](StageContext& context, const Json& config) {
        const AdaptiveLambdaOptions options = parse_options(config);
        const auto initial = exact_audit(context.db, context.density);
        const auto stats = adaptive::run_adaptive_lambda_search(
            context.db, context.density.bins_x, context.density.bins_y,
            context.density.target_density, options, options.verbose);
        if (std::abs(initial.hpwl - stats.initial.hpwl) >
                1.0e-6 * std::max(1.0, initial.hpwl) ||
            std::abs(initial.overflow_ratio - stats.initial.overflow) >
                1.0e-12) {
            throw std::runtime_error(
                "adaptive_lambda_gp internal audit disagrees with the stage-boundary audit");
        }
        print_metrics("initial", stats.initial);
        print_metrics("best_any", stats.best_any);
        print_metrics("best_feasible", stats.best_feasible);
        print_metrics("last_before_restore", stats.last_before_restore);
        print_metrics("final_selected", stats.final_selected);
        std::cout << "[adaptive_lambda_gp] max_excursion_percent=" << std::setprecision(14)
                  << stats.maximum_observed_overflow * 100.0
                  << " iterations_above_10pct=" << stats.iterations_above_10_percent
                  << " first_return_iter=" << stats.first_return_to_feasible_iteration
                  << " optimizer_resets=" << stats.optimizer_resets
                  << " lambda_initial=" << stats.lambda_initial
                  << " lambda_max=" << stats.lambda_max_observed
                  << " lambda_final=" << stats.lambda_final << '\n';
        return StageStats{options.iterations, stats.accepted, stats.rejected,
                          stats.objective_evaluations};
    });
}

}  // namespace nsgp::modules
