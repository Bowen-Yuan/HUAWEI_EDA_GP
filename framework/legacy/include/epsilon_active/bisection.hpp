#pragma once

#include "epsilon_active/density.hpp"
#include "epsilon_active/types.hpp"

namespace ea {

struct BisectionConfig {
    bool enabled = false;
    int leaf_bins = 8;
    int max_depth = -1;
    int degree_limit = 64;
    int fm_passes = 2;
    Real balance_tolerance = 0.02;
    bool position_seeded = false;
    // Preserve the current side of a cut whenever possible; only move
    // boundary surplus needed to satisfy the two exact capacity upper bounds.
    bool surplus_only = false;
    bool leaf_curve_order = false;
    bool leaf_nearest_capacity = false;
    bool leaf_hpwl_guided = false;
    // Optional locality guard for HPWL-guided leaf assignment.  A negative
    // value preserves the original unrestricted exact-capacity scan; a
    // nonnegative value limits candidate bins to a Chebyshev radius around
    // the node's current bin, preserving net geometry during continuation.
    int leaf_local_radius = -1;
    bool atomic_coarsening = false;
    int atomic_max_nodes = 8;
    int atomic_degree_limit = 16;
    int atomic_rounds = 3;
    int atomic_release_depth = 12;
    bool axis_rank_map = false;
    bool axis_map_only = false;
    bool curve_rank_map = false;
    bool curve_map_only = false;
};

struct BisectionStats {
    int splits = 0;
    int leaves = 0;
    int max_depth = 0;
    int fm_moves = 0;
    int atomic_groups = 0;
    int atomic_merges = 0;
    int atomic_unit_moves = 0;
    Real cut_net_weight = 0.0;
    Real maximum_capacity_error = 0.0;
    Real initial_hpwl = 0.0;
    Real initial_overflow = 0.0;
    Real final_hpwl = 0.0;
    Real final_overflow = 0.0;
    double wall_seconds = 0.0;
};

BisectionStats recursive_hypergraph_bisection(
    Database& db, ExactOverlapDensity& density, const BisectionConfig& config);

}  // namespace ea
