#include "microkernel.hpp"
#include "epsilon_active/placer.hpp"
#include "epsilon_active/lambda_controller.hpp"
#include "epsilon_active/optimizer.hpp"
#include "epsilon_active/epsilon_continuation.hpp"
#include <algorithm>
#include <string>

namespace nsgp::modules {
namespace {
double nested(const Json& c, const char* object, const char* key, double fallback) {
    return c.contains(object) && c.at(object).contains(key)
        ? c.at(object).at(key).get<double>() : fallback;
}
int nested_int(const Json& c, const char* object, const char* key, int fallback) {
    return c.contains(object) && c.at(object).contains(key)
        ? c.at(object).at(key).get<int>() : fallback;
}
}

// DREAMPlace-inspired overflow-driven epsilon continuation: the direction
// oracle radius widens under high overflow and tightens back near the
// target.  Exact metrics are untouched; optimizer resets only on band
// changes (state restart, never a layout rollback).
StageStats run_epsilon_continuation(StageContext& context, const Json& c) {
    ea::PlaceConfig p;
    p.bins_x = context.density.bins_x;
    p.bins_y = context.density.bins_y;
    p.target_density = context.density.target_density;
    p.threads = context.threads;
    p.iterations = c.value("iterations", 100);
    p.active_power = nested(c, "hpwl_direction", "active_power", 4.0);
    p.degree_limit = nested_int(c, "hpwl_direction", "degree_limit", 32);
    if (c.contains("hpwl_direction") &&
        c.at("hpwl_direction").contains("high_degree_mode")) {
        p.high_degree_mode = c.at("hpwl_direction").at("high_degree_mode").get<std::string>();
    }
    p.preconditioner = c.value("preconditioner", p.preconditioner);
    if (c.contains("optimizer")) {
        if (c.at("optimizer").contains("name")) {
            p.optimizer = ea::parse_optimizer(
                c.at("optimizer").at("name").get<std::string>());
        }
        p.step_fraction = nested(c, "optimizer", "learning_rate", 0.002);
        p.max_step_multiplier = nested(c, "optimizer", "maximum_delta", 0.25);
        p.beta1 = nested(c, "optimizer", "beta1", p.beta1);
        p.beta2 = nested(c, "optimizer", "beta2", p.beta2);
    }
    if (c.contains("epsilon_continuation")) {
        const auto& e = c.at("epsilon_continuation");
        p.epsilon_continuation = true;
        p.epsilon_schedule.overflow_low = e.value("overflow_low", 0.07);
        p.epsilon_schedule.overflow_high = e.value("overflow_high", 0.15);
        p.epsilon_schedule.wire_min_bins = e.value("wire_min_bins", 0.25);
        p.epsilon_schedule.wire_max_bins = e.value("wire_max_bins", 1.0);
        p.epsilon_schedule.density_min_bins = e.value("density_min_bins", 0.0);
        p.epsilon_schedule.density_max_bins = e.value("density_max_bins", 0.5);
        p.epsilon_schedule.restart_on_band_change =
            e.value("restart_on_band_change", true);
    }
    p.lambda.density_weight_scale = nested(c, "lambda", "density_weight_scale", 1.0);
    if (c.contains("lambda") && c.at("lambda").contains("policy")) {
        p.lambda.policy = ea::parse_lambda_policy(
            c.at("lambda").at("policy").get<std::string>());
    }
    p.lambda.update_interval = nested_int(c, "lambda", "update_interval", p.lambda.update_interval);
    p.lambda.trajectory_horizon = nested_int(c, "lambda", "trajectory_horizon", p.iterations);
    p.lambda.stop_overflow = nested(c, "lambda", "stop_overflow", 0.07);
    p.stop_overflow = nested(c, "selector", "overflow_cap", 0.07);
    p.log_every = c.value("log_every", 1);
    p.adaptive_epsilon = false;
    if (context.log_iteration) {
        p.iteration_hook = [&context](int iteration, ea::Real hpwl,
                                      ea::Real overflow_ratio, ea::Real lambda,
                                      ea::Real learning_rate) {
            context.log_iteration(iteration, "overflow_epsilon_continuation_gp", Json{
                {"hpwl", hpwl},
                {"overflow_percent", overflow_ratio * 100.0},
                {"lambda", lambda},
                {"learning_rate", learning_rate},
                {"objective_evaluations", 2}});
        };
    }
    const auto result = ea::global_place(context.db, p);
    clamp_movable(context.db);
    return StageStats{p.iterations, result.accepted_batches,
                      result.rejected_batches, result.objective_evaluations, 0, 0};
}
void register_overflow_epsilon_continuation_gp(nsgp::ModuleRegistry& registry) {
    registry.add("overflow_epsilon_continuation_gp", run_epsilon_continuation);
}
}  // namespace nsgp::modules
