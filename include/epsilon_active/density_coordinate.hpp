#pragma once

#include "epsilon_active/density.hpp"
#include "epsilon_active/types.hpp"

#include <vector>

namespace ea {

// Exact-overlap coordinate descent used as an optional active-set escape
// phase.  It never changes the HPWL/overlap oracles: a move is accepted only
// after exact rectangle-to-bin evaluation, with an optional cumulative HPWL
// budget acting as a search guard.
struct DensityCoordinateConfig {
    int sweeps = 0;
    int line_search_steps = 4;
    Real step_bins = 4.0;
    bool breakpoint_oracle = true;
    bool axis_separated = true;
    bool node_moves = true;
    Real hpwl_budget_fraction = 0.0;
    Real overflow_band_fraction = 0.0;
    int degree_limit = 100;
    bool net_blocks = false;
    bool capacity_assignment = false;
    int assignment_bins = 64;
    bool hierarchical_clusters = false;
    // Allow each node in a connected cluster to choose an exact breakpoint
    // fraction independently; this is an operator change, not smoothing.
    bool elastic_clusters = false;
    bool cooperative_batch = false;
    bool heavy_edge_growth = false;
    // Restrict hierarchical cluster destinations to 4-neighbour coarse bins.
    // This is an exact candidate-generation restriction, not a density
    // surrogate or smoothing term.
    bool adjacent_bin_only = false;
    int adjacent_bin_radius = 1;
    int hierarchy_levels = 3;
    int cluster_max_nodes = 128;
    int cluster_candidate_bins = 8;
    int cluster_max_clusters = 2000;
    int block_degree_limit = 16;
    int block_max_nodes = 32;
    int block_max_blocks = 20000;
};

struct DensityCoordinateSweepStats {
    int sweep = 0;
    int moves = 0;
    Real hpwl = 0.0;
    Real overflow = 0.0;
};

struct DensityCoordinateStats {
    int sweeps = 0;
    int moves = 0;
    int candidates = 0;
    int cluster_proposals = 0;
    int cluster_moves = 0;
    int cooperative_batch_trials = 0;
    int cooperative_batches = 0;
    Real initial_hpwl = 0.0;
    Real initial_overflow = 0.0;
    Real final_hpwl = 0.0;
    Real final_overflow = 0.0;
    double wall_seconds = 0.0;
    std::vector<DensityCoordinateSweepStats> trajectory;
};

DensityCoordinateStats coordinate_descent_overlap(
    Database& db, ExactOverlapDensity& density,
    const DensityCoordinateConfig& config);

}  // namespace ea
