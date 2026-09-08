#include "nonlocal_common.hpp"
#include "epsilon_active/hpwl.hpp"
#include <algorithm>
#include <cmath>

namespace nsgp::modules {
void register_finite_radius_oracle_gp(ModuleRegistry& r) {
    r.add("finite_radius_oracle_gp", [](StageContext& x,const Json& j) {
        const int probes=j.value("probe_count",32); const double scale=j.value("radius_bin_scale",1.0);
        return nonlocal::run(x,j,[probes,scale](const ea::Database& db0,ea::ExactOverlapDensity& d,const ea::AuxiliaryContext& ctx,
                                                std::vector<ea::Real>& ax,std::vector<ea::Real>& ay) {
            // The callback receives a const database by contract.  A one-sided exact
            // probe is performed on a short-lived copy, never on the live layout.
            ax.assign(db0.nodes.size(),0);ay.assign(db0.nodes.size(),0);
            const int count=std::min<int>(probes,db0.movable_ids.size()); const auto delta=scale*std::min(d.bin_width(),d.bin_height());
            const auto base_d=d.evaluate(0,1,nullptr,nullptr);
            std::vector<int> ids=db0.movable_ids;
            std::sort(ids.begin(),ids.end(),[&](int a,int b){const auto ra=ctx.wire_x[a]*ctx.wire_x[a]+ctx.wire_y[a]*ctx.wire_y[a]+ctx.exact_density_x[a]*ctx.exact_density_x[a]+ctx.exact_density_y[a]*ctx.exact_density_y[a];const auto rb=ctx.wire_x[b]*ctx.wire_x[b]+ctx.wire_y[b]*ctx.wire_y[b]+ctx.exact_density_x[b]*ctx.exact_density_x[b]+ctx.exact_density_y[b]*ctx.exact_density_y[b];return ra<rb;});
            for(int k=0;k<count;++k){const int id=ids[k];for(int axis=0;axis<2;++axis){double best=0,bestv=1e300;for(int sign:{-1,1}){ea::Database trial=db0;auto& n=trial.nodes[id];if(axis==0)n.x=std::clamp(n.x+sign*delta,trial.xl+.5*n.width,trial.xh-.5*n.width);else n.y=std::clamp(n.y+sign*delta,trial.yl+.5*n.height,trial.yh-.5*n.height);ea::ExactOverlapDensity den(trial,d.bins_x(),d.bins_y(),d.target_density());const auto value=den.evaluate(0,1,nullptr,nullptr).energy;if(value<bestv){bestv=value;best=sign*(value-base_d.energy)/delta;}}if(axis==0)ax[id]=best;else ay[id]=best;}}
        });
    });
}
}
