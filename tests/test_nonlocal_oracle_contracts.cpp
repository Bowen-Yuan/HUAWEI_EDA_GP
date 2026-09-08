#include "epsilon_active/nonlocal_descent.hpp"

#include <cmath>
#include <iostream>

int main() {
    ea::Database db;
    db.xl=0; db.yl=0; db.xh=20; db.yh=20; db.movable_area=16;
    db.nodes={{0,"m",5,5,4,4,false,false,"N"}};
    db.movable_ids={0}; db.node_pin_count={0}; db.node_pin_offsets={0,0};
    ea::NonlocalDescentConfig c; c.iterations=1; c.learning_rate=.002;
    c.maximum_delta=.25; c.auxiliary_density_weight=.1;
    const auto result=ea::run_nonlocal_descent(db,2,2,1.0,c,
        [](const ea::Database& layout,ea::ExactOverlapDensity&,const ea::AuxiliaryContext&,
           std::vector<ea::Real>& x,std::vector<ea::Real>& y) {
            x.assign(layout.nodes.size(),0); y.assign(layout.nodes.size(),0);
            x[0]=-1.0; // Gradient negative means the descent update moves right.
        });
    // A zero wire/density baseline must remain finite and never trigger an
    // accept/reject path merely because the auxiliary is sparse.
    if (!std::isfinite(db.nodes[0].x) || !std::isfinite(db.nodes[0].y)) return 1;
    if (result.objective_evaluations < 2) return 2;
    std::cout << "nonlocal oracle contracts passed\n";
}
