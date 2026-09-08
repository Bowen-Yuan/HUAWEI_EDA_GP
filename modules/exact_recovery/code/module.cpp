#include "microkernel.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/recovery.hpp"

namespace nsgp::modules {
void register_exact_recovery(ModuleRegistry& registry) {
    registry.add("exact_recovery", [](StageContext& context, const Json& c) {
        ea::ExactOverlapDensity d(context.db,context.density.bins_x,
                                  context.density.bins_y,context.density.target_density);
        ea::RecoveryConfig r; r.sweeps=c.value("sweeps",1);
        r.line_search_steps=c.value("line_search_steps",4);
        r.step_bins=c.value("step_bins",8.0); r.overflow_cap=c.value("overflow_cap",.07);
        r.degree_limit=c.value("degree_limit",100);
        r.axis_separated=c.value("axis_separated",false);
        r.breakpoint_oracle=c.value("breakpoint_oracle",false);
        r.net_block=c.value("net_block",false);
        r.net_block_density_direction=c.value("net_block_density_direction",false);
        r.net_block_degree_limit=c.value("net_block_degree_limit",16);
        r.net_block_max_nodes=c.value("net_block_max_nodes",8);
        r.net_block_max_blocks=c.value("net_block_max_blocks",10000);
        r.compact_directions=c.value("compact_directions",false);
        const auto result=ea::recover_hpwl_under_overflow(context.db,d,r);
        return StageStats{result.sweeps,result.moves,0,result.objective_evaluations};
    });
}
}
