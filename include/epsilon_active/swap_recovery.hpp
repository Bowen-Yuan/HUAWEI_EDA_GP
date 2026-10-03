#pragma once

#include "epsilon_active/density.hpp"
#include "epsilon_active/types.hpp"

#include <vector>

namespace ea {

struct SwapRecoveryConfig {
    int sweeps = 0;
    int radius_bins = 4;
    int candidates = 16;
    int exact_shortlist = 0;
    int degree_limit = 100;
    int permutation_sweeps = 0;
    int assignment_sweeps = 0;
    int assignment_radius_bins = 16;
    int assignment_candidates = 8;
    int assignment_max_bids_per_node = 32;
    Real assignment_epsilon = 1.0;
    // Optional density-guided equal-shape swaps.  A positive HPWL budget is
    // fractional relative to the sweep-start exact HPWL and is paid only
    // when exact overflow decreases.
    bool density_guided = false;
    Real density_hpwl_budget = 0.0;
    Real density_score_weight = 1.0;
    bool allow_unequal_shapes = false;
    // Prioritize exact equal-shape exchanges with nodes sharing a net before
    // the spatial candidate ring. This changes only the search neighborhood;
    // occupancy and the audited nonsmooth objective are unchanged.
    bool net_aware = false;
};

struct SwapSweepStats {
    int sweep = 0;
    int swaps = 0;
    Real hpwl = 0.0;
    Real overflow = 0.0;
};

struct SwapRecoveryStats {
    int sweeps = 0;
    int swaps = 0;
    int permutation_sweeps = 0;
    int permuted_nodes = 0;
    int assignment_sweeps = 0;
    int assignment_batches = 0;
    int assigned_nodes = 0;
    std::int64_t assignment_bids = 0;
    int objective_evaluations = 0;
    Real initial_hpwl = 0.0;
    Real initial_overflow = 0.0;
    Real final_hpwl = 0.0;
    Real final_overflow = 0.0;
    double wall_seconds = 0.0;
    std::vector<SwapSweepStats> trajectory;
};

SwapRecoveryStats recover_hpwl_with_equal_shape_swaps(
    Database& db, ExactOverlapDensity& density,
    const SwapRecoveryConfig& config);

}  // namespace ea
