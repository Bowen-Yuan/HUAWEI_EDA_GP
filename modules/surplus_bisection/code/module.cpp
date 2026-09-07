#include "microkernel.hpp"
#include "epsilon_active/bisection.hpp"
#include "epsilon_active/density.hpp"

namespace nsgp::modules {
void register_surplus_bisection(ModuleRegistry& registry) {
    registry.add("surplus_bisection", [](StageContext& context, const Json& c) {
        ea::ExactOverlapDensity d(context.db,context.density.bins_x,
                                  context.density.bins_y,context.density.target_density);
        ea::BisectionConfig b; b.enabled=true; b.leaf_bins=c.value("leaf_bins",8);
        b.position_seeded=c.value("position_seeded",true);
        b.surplus_only=c.value("surplus_only",true);
        b.leaf_nearest_capacity=c.value("leaf_nearest_capacity",true);
        b.leaf_hpwl_guided=c.value("leaf_hpwl_guided",true);
        const auto result=ea::recursive_hypergraph_bisection(context.db,d,b);
        return StageStats{1,result.fm_moves,0,4};
    });
}
}
