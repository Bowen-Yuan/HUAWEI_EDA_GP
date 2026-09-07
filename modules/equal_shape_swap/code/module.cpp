#include "microkernel.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/swap_recovery.hpp"

namespace nsgp::modules {
void register_equal_shape_swap(ModuleRegistry& registry) {
    registry.add("equal_shape_swap", [](StageContext& context, const Json& c) {
        ea::ExactOverlapDensity d(context.db,context.density.bins_x,
                                  context.density.bins_y,context.density.target_density);
        ea::SwapRecoveryConfig s; s.sweeps=c.value("sweeps",1);
        s.radius_bins=c.value("radius_bins",4); s.candidates=c.value("candidates",16);
        s.exact_shortlist=c.value("exact_shortlist",0);
        s.net_aware=c.value("net_aware",true);
        const auto result=ea::recover_hpwl_with_equal_shape_swaps(context.db,d,s);
        return StageStats{result.sweeps,result.swaps,0,result.objective_evaluations};
    });
}
}
