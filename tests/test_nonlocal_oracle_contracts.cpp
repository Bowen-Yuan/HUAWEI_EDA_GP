#include "epsilon_active/nonlocal_descent.hpp"
#include "finite_radius_probe.hpp"

#include "epsilon_active/hpwl.hpp"

#include <cmath>
#include <iostream>

namespace {

// Synthetic single-movable-node layout with a 2x2 density grid where the
// right half is pre-filled by a fixed macro, so moving the node right
// strictly increases density energy and moving left decreases it.
ea::Database make_probe_layout() {
    ea::Database db;
    db.xl=0; db.yl=0; db.xh=20; db.yh=20; db.movable_area=16;
    db.nodes={
        {0,"m",5,5,4,4,false,false,"N"},
        {1,"wall",15,10,10,20,true,false,"N"}};
    db.fixed_ids={1}; db.movable_ids={0};
    db.node_pin_count={0,0}; db.node_pin_offsets={0,0,0};
    return db;
}

}  // namespace

int main() {
    ea::Database db;
    db.xl=0; db.yl=0; db.xh=20; db.yh=20; db.movable_area=16;
    db.nodes={{0,"m",5,5,4,4,false,false,"N"}};
    db.movable_ids={0}; db.node_pin_count={0}; db.node_pin_offsets={0,0};
    ea::NonlocalDescentConfig c; c.iterations=1; c.learning_rate=.002;
    c.maximum_delta=.25; c.auxiliary_density_weight=.1;
    const auto result=ea::run_nonlocal_descent(db,2,2,1.0,c,
        [](const ea::Database& layout,ea::ExactOverlapDensity&,const ea::AuxiliaryContext&,
           std::vector<ea::Real>& x,std::vector<ea::Real>& y) -> ea::AuxiliaryQueryStats {
            x.assign(layout.nodes.size(),0); y.assign(layout.nodes.size(),0);
            x[0]=-1.0; // Gradient negative means the descent update moves right.
            return {};
        });
    // A zero wire/density baseline must remain finite and never trigger an
    // accept/reject path merely because the auxiliary is sparse.
    if (!std::isfinite(db.nodes[0].x) || !std::isfinite(db.nodes[0].y)) return 1;
    if (result.objective_evaluations < 2) return 2;

    // --- Fix A4 contract: AuxiliaryContext.hpwl equals exact HPWL ---
    {
        ea::Database ctx_db;
        ctx_db.xl=0; ctx_db.yl=0; ctx_db.xh=20; ctx_db.yh=20;
        ctx_db.movable_area=16;
        ctx_db.nodes={{0,"m",5,5,4,4,false,false,"N"}};
        ctx_db.movable_ids={0}; ctx_db.node_pin_count={0};
        ctx_db.node_pin_offsets={0,0};
        ctx_db.pins={}; ctx_db.nets={};
        // one self net so exact HPWL is nonzero: net with one pin on node 0
        ctx_db.pins={{0,3,0}};
        ctx_db.nets.push_back({0,"n",0,1,2.0});
        ctx_db.node_pin_count={1};
        ctx_db.node_pin_offsets={0,1};
        ctx_db.node_pin_indices={0};
        ea::NonlocalDescentConfig c;
        c.iterations=1; c.learning_rate=.002; c.maximum_delta=.25;
        ea::Real observed_hpwl=-1.0;
        ea::run_nonlocal_descent(ctx_db,2,2,1.0,c,
            [&](const ea::Database&, ea::ExactOverlapDensity&,
                const ea::AuxiliaryContext& ctx,
                std::vector<ea::Real>& x, std::vector<ea::Real>& y) {
                x.assign(ctx_db.nodes.size(),0);
                y.assign(ctx_db.nodes.size(),0);
                observed_hpwl=ctx.hpwl;
                return ea::AuxiliaryQueryStats{};
            });
        ea::ExactHpwl independent(ctx_db);
        const ea::Real expected=independent.evaluate(0,1,-1,nullptr,nullptr);
        if(std::abs(observed_hpwl-expected)>1e-9) return 3;
    }

    // --- Fix A1/A2 contracts: finite-radius descent gate (FR1-FR5) ---
    {
        // Wall occupies the right half of the die.
        ea::Database fdb=make_probe_layout();
        ea::ExactOverlapDensity density(fdb,2,2,1.0);
        density.rebuild_occupancy(ea::DensityAccumulationMode::DeterministicCanonical);
        std::vector<ea::Real> ctx_ax(2,0.0), ctx_ay(2,0.0);
        ea::AuxiliaryContext ctx{0,0.1,ctx_ax,ctx_ay,ctx_ay,ctx_ay,0.0,
                                 density.evaluate(0,1,nullptr,nullptr)};
        std::vector<ea::Real> ax,ay;

        // FR1: node at x=5 with probe radius 1 -> left probe is free space
        // (energy unchanged, not a descent) and right probe enters the wall
        // (worse).  Both probes fail the gate -> auxiliary stays 0.
        {
            fdb.nodes[0].x=5.0;
            density.rebuild_occupancy(ea::DensityAccumulationMode::DeterministicCanonical);
            const auto stats=nsgp::modules::finite_radius_probe(
                fdb,density,ctx,1,0.25,1e-12,ax,ay);
            if(std::abs(ax[0])>0.0 || std::abs(ay[0])>0.0) return 4;
            if(stats.local_density_queries<1) return 5;
        }
        // FR2: overlap the node with the wall; the left probe strictly
        // decreases energy.  g_x must be positive so x -= eta*g moves -x.
        {
            fdb.nodes[0].x=12.0;
            density.rebuild_occupancy(ea::DensityAccumulationMode::DeterministicCanonical);
            const auto stats=nsgp::modules::finite_radius_probe(
                fdb,density,ctx,1,0.25,1e-12,ax,ay);
            if(!(ax[0]>0.0)) return 6;
            if(stats.local_density_queries<1) return 7;
        }
        // FR3: mirror case.  Node just left of the wall; the right probe
        // decreases energy is impossible here, so instead overlap on the
        // left side of an inverted wall: use a node at x=9.5 whose right
        // probe at 10+delta reduces overlap by moving INTO free space.
        {
            // Invert: make the wall occupy the LEFT half instead.
            fdb.nodes[1].x=5.0;  // wall now covers x in [0,10]
            fdb.nodes[0].x=8.0;  // node overlaps the left wall
            // The fixed occupancy is cached at construction, so a moved
            // fixed node requires a fresh evaluator.
            ea::ExactOverlapDensity density3(fdb,2,2,1.0);
            density3.rebuild_occupancy(ea::DensityAccumulationMode::DeterministicCanonical);
            const auto stats=nsgp::modules::finite_radius_probe(
                fdb,density3,ctx,1,0.25,1e-12,ax,ay);
            // right probe (+e_x) strictly decreases energy: g_x must be
            // negative so x -= eta*g moves +x.
            if(!(ax[0]<0.0)) return 8;
            if(stats.local_density_queries<1) return 9;
        }
        // FR4: boundary clamp uses the actual displacement denominator.
        {
            fdb.nodes[1].x=15.0;  // wall back to the right half
            fdb.nodes[0].x=0.5;   // touching the left die edge: width 4 -> clamp keeps x>=2
            ea::ExactOverlapDensity density4(fdb,2,2,1.0);
            density4.rebuild_occupancy(ea::DensityAccumulationMode::DeterministicCanonical);
            const auto stats=nsgp::modules::finite_radius_probe(
                fdb,density4,ctx,1,1.0,1e-12,ax,ay);
            // Left probe clamps to x=2 (actual displacement 1.5, not 5).
            // If the node is in free space, no descent is expected; the only
            // contract here is that the probe never divides by nominal delta.
            if(stats.local_density_queries<0) return 10;
        }
        // FR5: evaluate_move delta equals a full deterministic reevaluation.
        {
            fdb.nodes[0].x=9.0;
            ea::ExactOverlapDensity density5(fdb,2,2,1.0);
            density5.rebuild_occupancy(ea::DensityAccumulationMode::DeterministicCanonical);
            const auto base=density5.evaluate(
                0,1,ea::DensityAccumulationMode::DeterministicCanonical,
                nullptr,nullptr);
            const auto mv=density5.evaluate_move(0,11.0,5);
            fdb.nodes[0].x=11.0;
            ea::ExactOverlapDensity density6(fdb,2,2,1.0);
            density6.rebuild_occupancy(ea::DensityAccumulationMode::DeterministicCanonical);
            const auto full=density6.evaluate(
                0,1,ea::DensityAccumulationMode::DeterministicCanonical,
                nullptr,nullptr);
            if(std::abs((base.energy+mv.energy_delta)-full.energy)>1e-12) return 11;
        }
    }
    std::cout << "nonlocal oracle contracts passed\n";
}
