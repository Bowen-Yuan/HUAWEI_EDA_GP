#include "microkernel.hpp"
#include "epsilon_active/placer.hpp"
#include <algorithm>

namespace nsgp::modules {
namespace {
double nested(const Json& c,const char* object,const char* key,double fallback) {
    return c.contains(object) && c.at(object).contains(key)
        ? c.at(object).at(key).get<double>() : fallback;
}
}
void register_exact_joint_gp(ModuleRegistry& registry) {
    registry.add("exact_joint_gp", [](StageContext& context, const Json& c) {
        ea::PlaceConfig p;
        p.bins_x=context.density.bins_x; p.bins_y=context.density.bins_y;
        p.target_density=context.density.target_density; p.threads=context.threads;
        p.iterations=c.value("iterations",10);
        p.hpwl_epsilon=nested(c,"hpwl_direction","epsilon",c.value("epsilon",125.0));
        p.active_power=nested(c,"hpwl_direction","active_power",c.value("active_power",4.0));
        p.degree_limit=c.contains("hpwl_direction")
            ? c.at("hpwl_direction").value("degree_limit",c.value("degree_limit",100))
            : c.value("degree_limit",100);
        p.step_fraction=nested(c,"optimizer","learning_rate",c.value("learning_rate",.003));
        p.max_step_multiplier=nested(c,"optimizer","maximum_delta",c.value("maximum_delta",4.0));
        p.lambda.density_weight_scale=nested(c,"lambda","density_weight_scale",.5);
        p.stop_overflow=nested(c,"selector","overflow_cap",c.value("overflow_cap",.07));
        if (c.contains("exact_batch_acceptance")) {
            p.batch_acceptance.enabled=c.at("exact_batch_acceptance").value("enabled",false);
            p.batch_acceptance.max_trials=c.at("exact_batch_acceptance").value("max_backtracks",0);
        }
        p.snapshot_every=0; p.output_dir.clear();
        p.coarse_flow.bins_x=std::min(64,context.density.bins_x);
        p.coarse_flow.bins_y=std::min(64,context.density.bins_y);
        const auto result=ea::global_place(context.db,p); clamp_movable(context.db);
        return StageStats{p.iterations,result.accepted_batches,result.rejected_batches,
                          result.objective_evaluations};
    });
}
}
