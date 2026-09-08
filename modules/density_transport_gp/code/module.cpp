#include "nonlocal_common.hpp"
#include <algorithm>
#include <limits>

namespace nsgp::modules {
void register_density_transport_gp(ModuleRegistry& r) {
    r.add("density_transport_gp", [](StageContext& x,const Json& j) {
        const int coarse=j.value("coarse_bins",32);
        return nonlocal::run(x,j,[coarse](const ea::Database& db,ea::ExactOverlapDensity& d,const ea::AuxiliaryContext&,
                                          std::vector<ea::Real>& ax,std::vector<ea::Real>& ay) {
            const int nx=d.bins_x(),ny=d.bins_y(),cx=std::min(coarse,nx),cy=std::min(coarse,ny);
            std::vector<ea::Real> mass(cx*cy,0),capacity(cx*cy,0); const auto& occ=d.occupancy();
            for(int y=0;y<ny;++y)for(int xx=0;xx<nx;++xx){int X=xx*cx/nx,Y=y*cy/ny;mass[Y*cx+X]+=occ[y*nx+xx];capacity[Y*cx+X]+=d.target_density()*d.bin_area();}
            std::vector<int> sink(cx*cy,-1);
            for(int s=0;s<cx*cy;++s)if(mass[s]>capacity[s]){ea::Real best=std::numeric_limits<ea::Real>::max();for(int t=0;t<cx*cy;++t)if(mass[t]<capacity[t]){const auto cost=std::abs(s%cx-t%cx)+std::abs(s/cx-t/cx);if(cost<best){best=cost;sink[s]=t;}}}
            ax.assign(db.nodes.size(),0);ay.assign(db.nodes.size(),0);
            for(int id:db.movable_ids){const auto& n=db.nodes[id];int X=std::clamp(int((n.x-db.xl)/(db.xh-db.xl)*cx),0,cx-1),Y=std::clamp(int((n.y-db.yl)/(db.yh-db.yl)*cy),0,cy-1);int s=Y*cx+X;if(sink[s]<0)continue;int t=sink[s];ax[id]=double(X-t%cx);ay[id]=double(Y-t/cx);}
        });
    });
}
}
