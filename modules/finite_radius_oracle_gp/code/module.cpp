#include "nonlocal_common.hpp"
#include "finite_radius_probe.hpp"

namespace nsgp::modules {
void register_finite_radius_oracle_gp(ModuleRegistry& r) {
    r.add("finite_radius_oracle_gp", [](StageContext& x,const Json& j) {
        const int probes=j.value("probe_count",32); const double scale=j.value("radius_bin_scale",1.0);
        const double descent_tolerance=j.value("descent_tolerance",1e-12);
        if (probes < 0 || scale <= 0.0 || descent_tolerance < 0.0)
            throw std::invalid_argument("finite_radius_oracle_gp requires probe_count >= 0, radius_bin_scale > 0, descent_tolerance >= 0");
        return nonlocal::run(x,j,"finite_radius_oracle_gp",
            [probes,scale,descent_tolerance](const ea::Database& db0,ea::ExactOverlapDensity& d,
                                             const ea::AuxiliaryContext& ctx,
                                             std::vector<ea::Real>& ax,std::vector<ea::Real>& ay) {
                return finite_radius_probe(db0,d,ctx,probes,scale,descent_tolerance,ax,ay);
            });
    });
}
}
