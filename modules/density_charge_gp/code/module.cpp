#include "nonlocal_common.hpp"
#include <algorithm>
#include <cmath>

namespace nsgp::modules {
void register_density_charge_gp(ModuleRegistry& r) {
    r.add("density_charge_gp", [](StageContext& x,const Json& j) {
        const double radius_scale=j.value("radius_bin_scale",4.0);
        return nonlocal::run(x,j,[radius_scale](const ea::Database& db,ea::ExactOverlapDensity& d,const ea::AuxiliaryContext&,
                                                std::vector<ea::Real>& ax,std::vector<ea::Real>& ay) {
            const int nx=d.bins_x(),ny=d.bins_y(); const auto& occ=d.occupancy();const auto cap=d.target_density()*d.bin_area();
            const auto h=std::min(d.bin_width(),d.bin_height());const auto radius=radius_scale*h;
            ax.assign(db.nodes.size(),0);ay.assign(db.nodes.size(),0);
            for(int id:db.movable_ids){const auto& n=db.nodes[id];int cx=std::clamp(int((n.x-db.xl)/d.bin_width()),0,nx-1),cy=std::clamp(int((n.y-db.yl)/d.bin_height()),0,ny-1);
                const int rx=std::max(1,int(radius/d.bin_width())+1),ry=std::max(1,int(radius/d.bin_height())+1);
                for(int yy=std::max(0,cy-ry);yy<=std::min(ny-1,cy+ry);++yy)for(int xx=std::max(0,cx-rx);xx<=std::min(nx-1,cx+rx);++xx){
                    const auto bx=db.xl+(xx+.5)*d.bin_width(),by=db.yl+(yy+.5)*d.bin_height();const auto dx=n.x-bx,dy=n.y-by;const auto dist=std::abs(dx)+std::abs(dy);if(dist>=radius)continue;
                    const auto q=occ[yy*nx+xx]-cap; ax[id]-=q*(dx>=0?1:-1);ay[id]-=q*(dy>=0?1:-1); }
            }
        });
    });
}
}
