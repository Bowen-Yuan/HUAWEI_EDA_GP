#include "nonlocal_common.hpp"
#include "epsilon_active/hpwl.hpp"
#include <algorithm>
#include <cmath>

namespace nsgp::modules {
void register_finite_radius_oracle_gp(ModuleRegistry& r) {
    r.add("finite_radius_oracle_gp", [](StageContext& x,const Json& j) {
        const int probes=j.value("probe_count",32); const double scale=j.value("radius_bin_scale",1.0);
        return nonlocal::run(x,j,[probes,scale](const ea::Database& db0,ea::ExactOverlapDensity& d,
                                                std::vector<ea::Real>& ax,std::vector<ea::Real>& ay) {
            // The callback receives a const database by contract.  A one-sided exact
            // probe is performed on a short-lived copy, never on the live layout.
            ax.assign(db0.nodes.size(),0);ay.assign(db0.nodes.size(),0);
            const int count=std::min<int>(probes,db0.movable_ids.size()); const auto delta=scale*std::min(d.bin_width(),d.bin_height());
            ea::ExactHpwl base_hpwl(db0);const auto base_d=d.evaluate(0,1,nullptr,nullptr);const auto base_w=base_hpwl.evaluate(0,1,-1,nullptr,nullptr);
            for(int k=0;k<count;++k){const int id=db0.movable_ids[k];for(int axis=0;axis<2;++axis){ea::Database trial=db0;auto& n=trial.nodes[id];if(axis==0)n.x=std::clamp(n.x+delta,trial.xl+.5*n.width,trial.xh-.5*n.width);else n.y=std::clamp(n.y+delta,trial.yl+.5*n.height,trial.yh-.5*n.height);ea::ExactHpwl h(trial);ea::ExactOverlapDensity den(trial,d.bins_x(),d.bins_y(),d.target_density());const auto dm=den.evaluate(0,1,nullptr,nullptr);const auto w=h.evaluate(0,1,-1,nullptr,nullptr);const auto secant=((w-base_w)/std::max(1.0,std::abs(base_w))+(dm.energy-base_d.energy)/std::max(1.0,std::abs(base_d.energy)+1.0))/delta;if(axis==0)ax[id]=secant;else ay[id]=secant;}}
        });
    });
}
}
