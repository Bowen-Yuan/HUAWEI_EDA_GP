#include "microkernel.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/transport.hpp"

namespace nsgp::modules {
void register_global_capacity_transport(ModuleRegistry& registry) {
    registry.add("global_capacity_transport", [](StageContext& context, const Json& c) {
        ea::ExactOverlapDensity d(context.db,context.density.bins_x,
                                  context.density.bins_y,context.density.target_density);
        ea::TransportConfig t;
        t.rounds=c.value("rounds",1); t.max_moves=c.value("max_moves",250000);
        t.max_source_bins=c.value("max_source_bins",64);
        t.candidate_lookahead=c.value("candidate_lookahead",8);
        t.group_size=c.value("group_size",1);
        t.destination_mode=ea::parse_transport_destination_mode(
            c.value("destination_mode",std::string("nearest")));
        t.global_hpwl_budget=c.value("global_hpwl_budget",-1.0);
        const auto result=ea::transport_excess_to_capacity(context.db,d,t);
        return StageStats{result.rounds,result.moves,0,4};
    });
}
}
