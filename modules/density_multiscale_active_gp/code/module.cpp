#include "nonlocal_common.hpp"
#include <cmath>

namespace nsgp::modules {
void register_density_multiscale_active_gp(ModuleRegistry& r) {
    r.add("density_multiscale_active_gp", [](StageContext& x,const Json& j) {
        const std::vector<ea::Real> scales=j.value("epsilon_bin_scales",std::vector<ea::Real>{.25,.5,1.,2.});
        return nonlocal::run(x,j,[scales](const ea::Database&,ea::ExactOverlapDensity& d,const ea::AuxiliaryContext&,
                                          std::vector<ea::Real>& ax,std::vector<ea::Real>& ay) {
            for (const auto s:scales) { std::vector<ea::Real> gx,gy;
                d.evaluate(s*std::min(d.bin_width(),d.bin_height()),1.,&gx,&gy);
                if (ax.size()!=gx.size()) { ax.assign(gx.size(),0); ay.assign(gy.size(),0); }
                double sum=0;for(std::size_t i=0;i<gx.size();++i)sum+=gx[i]*gx[i]+gy[i]*gy[i];const double r=std::sqrt(sum/std::max<std::size_t>(1,2*gx.size()));if(r<=1e-18)continue;
                for(std::size_t i=0;i<gx.size();++i){ax[i]+=gx[i]/r;ay[i]+=gy[i]/r;} }
        });
    });
}
}
