#include "microkernel.hpp"
#include "epsilon_active/placer.hpp"
#include <algorithm>

namespace nsgp::modules {
void register_hpwl_adam(ModuleRegistry& registry) {
    registry.add("hpwl_adam", [](StageContext& context, const Json& c) {
        ea::PlaceConfig p;
        p.bins_x=context.density.bins_x; p.bins_y=context.density.bins_y;
        p.target_density=context.density.target_density; p.threads=context.threads;
        p.iterations=c.value("iterations",10); p.hpwl_epsilon=c.value("epsilon",125.0);
        p.active_power=c.value("active_power",4.0); p.degree_limit=c.value("degree_limit",100);
        p.step_fraction=c.value("learning_rate",.003);
        p.max_step_multiplier=c.value("maximum_delta",4.0);
        p.lambda.density_weight_scale=0.0; p.snapshot_every=0; p.output_dir.clear();
        p.coarse_flow.bins_x=std::min(64,context.density.bins_x);
        p.coarse_flow.bins_y=std::min(64,context.density.bins_y);
        const auto result=ea::global_place(context.db,p); clamp_movable(context.db);
        return StageStats{p.iterations,result.accepted_batches,result.rejected_batches,
                          result.objective_evaluations};
    });
}
}
