#include "nonlocal_common.hpp"
#include <algorithm>

namespace nsgp::modules {
void register_density_cut_pressure_gp(ModuleRegistry& r) {
    r.add("density_cut_pressure_gp", [](StageContext& x,const Json& j) {
        return nonlocal::run(x,j,[](const ea::Database& db,ea::ExactOverlapDensity& d,const ea::AuxiliaryContext&,
                                    std::vector<ea::Real>& ax,std::vector<ea::Real>& ay) {
            const int nx=d.bins_x(), ny=d.bins_y(); const auto& occ=d.occupancy();
            std::vector<ea::Real> px(nx,0),py(ny,0); const ea::Real cap=d.target_density()*d.bin_area();
            for(int y=0;y<ny;++y)for(int xx=0;xx<nx;++xx){const auto q=occ[y*nx+xx]-cap;px[xx]+=q;py[y]+=q;}
            for(int i=1;i<nx;++i)px[i]+=px[i-1]; for(int i=1;i<ny;++i)py[i]+=py[i-1];
            ax.assign(db.nodes.size(),0); ay.assign(db.nodes.size(),0);
            for(int id:db.movable_ids){const auto& n=db.nodes[id];int bx=std::clamp(int((n.x-db.xl)/d.bin_width()),0,nx-1),by=std::clamp(int((n.y-db.yl)/d.bin_height()),0,ny-1);
                ax[id]=-(bx==0?0:px[bx-1]); ay[id]=-(by==0?0:py[by-1]);}
        });
    });
}
}
