#include "microkernel.hpp"
#include "epsilon_active/placer.hpp"
#include "epsilon_active/lambda_controller.hpp"
#include "epsilon_active/optimizer.hpp"
#include "epsilon_active/spectral_step.hpp"
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

// DREAMPlace-inspired spectral (Barzilai-Borwein) step-size stage.  The
// direction oracle stays the exact joint GP one; only the scalar learning
// rate is adapted from s/y curvature.  Pure-descent: no rollback.
StageStats run_spectral_bb(StageContext& context, const Json& c) {
    ea::PlaceConfig p;
    p.bins_x = context.density.bins_x;
    p.bins_y = context.density.bins_y;
    p.target_density = context.density.target_density;
    p.threads = context.threads;
    p.iterations = c.value("iterations", 100);
    p.hpwl_epsilon = nested(c, "hpwl_direction", "epsilon", 125.0);
    p.active_power = nested(c, "hpwl_direction", "active_power", 4.0);
    p.degree_limit = nested_int(c, "hpwl_direction", "degree_limit", 32);
    if (c.contains("hpwl_direction") &&
        c.at("hpwl_direction").contains("high_degree_mode")) {
        p.high_degree_mode = c.at("hpwl_direction").at("high_degree_mode").get<std::string>();
    }
    p.preconditioner = c.value("preconditioner", p.preconditioner);
    p.net_batch.weight = nested(c, "net_batch", "weight", 0.0);
    p.net_batch.degree_limit =
        nested_int(c, "net_batch", "degree_limit", 32);
    p.optimizer = ea::OptimizerKind::SGD;  // spectral owns the step scale
    p.spectral_step = true;
    if (c.contains("spectral_step")) {
        const auto& s = c.at("spectral_step");
        p.spectral.method = s.value("method", p.spectral.method);
        p.spectral.base_learning_rate = s.value("base_learning_rate", 0.002);
        p.spectral.alpha_min_ratio = s.value("alpha_min_ratio", p.spectral.alpha_min_ratio);
        p.spectral.alpha_max_ratio = s.value("alpha_max_ratio", p.spectral.alpha_max_ratio);
        p.spectral.max_step_growth = s.value("max_step_growth", p.spectral.max_step_growth);
        p.spectral.curvature_epsilon = s.value("curvature_epsilon", p.spectral.curvature_epsilon);
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
    p.max_step_multiplier = nested(c, "spectral_step", "maximum_delta", 0.25);
    if (context.log_iteration) {
        p.iteration_hook = [&context](int iteration, ea::Real hpwl,
                                      ea::Real overflow_ratio, ea::Real lambda,
                                      ea::Real learning_rate) {
            context.log_iteration(iteration, "spectral_bb_gp", Json{
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
void register_spectral_bb_gp(nsgp::ModuleRegistry& registry) {
    registry.add("spectral_bb_gp", run_spectral_bb);
}
}  // namespace nsgp::modules
