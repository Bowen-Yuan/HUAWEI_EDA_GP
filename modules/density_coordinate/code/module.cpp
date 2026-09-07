#include "microkernel.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/density_coordinate.hpp"

namespace nsgp::modules {
void register_density_coordinate(ModuleRegistry& registry) {
    registry.add("density_coordinate", [](StageContext& context, const Json& c) {
        ea::ExactOverlapDensity d(context.db,context.density.bins_x,
                                  context.density.bins_y,context.density.target_density);
        ea::DensityCoordinateConfig v;
        v.sweeps=c.value("sweeps",1);
        v.line_search_steps=c.value("line_search_steps",4);
        v.step_bins=c.value("step_bins",4.0);
        v.hpwl_budget_fraction=c.value("hpwl_budget_fraction",0.0);
        v.overflow_band_fraction=c.value("overflow_band_fraction",0.0);
        v.net_blocks=c.value("net_blocks",false);
        const auto result=ea::coordinate_descent_overlap(context.db,d,v);
        return StageStats{result.sweeps,result.moves,
                          result.candidates-result.moves,result.candidates};
    });
}
}
