#include "nonlocal_common.hpp"

namespace nsgp::modules {
void register_density_multiscale_active_gp(ModuleRegistry& r) {
    r.add("density_multiscale_active_gp", [](StageContext& x,const Json& j) {
        const std::vector<ea::Real> scales=j.value("epsilon_bin_scales",std::vector<ea::Real>{0,.25,.5,1.,2.});
        return nonlocal::run(x,j,[scales](const ea::Database&,ea::ExactOverlapDensity& d,
                                          std::vector<ea::Real>& ax,std::vector<ea::Real>& ay) {
            for (const auto s:scales) { std::vector<ea::Real> gx,gy;
                d.evaluate(s*std::min(d.bin_width(),d.bin_height()),1.,&gx,&gy);
                if (ax.size()!=gx.size()) { ax.assign(gx.size(),0); ay.assign(gy.size(),0); }
                for(std::size_t i=0;i<gx.size();++i){ax[i]+=gx[i]/scales.size();ay[i]+=gy[i]/scales.size();} }
        });
    });
}
}
