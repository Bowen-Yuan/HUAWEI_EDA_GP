#pragma once

#include "epsilon_active/density.hpp"
#include "epsilon_active/types.hpp"

#include <vector>
#include <cstdint>

namespace ea {

struct RecoveryConfig {
    int sweeps = 0;
    int line_search_steps = 4;
    int degree_limit = 100;
    Real step_bins = 8.0;
    Real hpwl_epsilon = 125.0;
    Real active_power = 4.0;
    Real overflow_cap = -1.0;
    // Optional exact bridge weight. Accepted score is HPWL delta plus this
    // weight times normalized exact overflow-area delta.
    Real objective_density_weight = 0.0;
    bool axis_separated = false;
    bool breakpoint_oracle = false;
    // Candidate-only active-set mechanisms.  They never alter exact overlap.
    bool boundary_active = false;
    Real boundary_epsilon_bins = 0.25;
    bool dead_zone_escape = false;
    bool controlled_exact_noise = false;
    std::uint64_t exact_noise_seed = 17;
    bool compact_directions = false;
    Real net_batch_weight = 0.0;
    int net_batch_degree_limit = 64;
    bool net_block = false;
    bool net_block_density_direction = false;
    bool net_block_contraction = false;
    bool node_moves = true;
    int net_block_degree_limit = 16;
    int net_block_max_nodes = 8;
    int net_block_max_blocks = 10000;
};

struct RecoverySweepStats {
    int sweep = 0;
    int moves = 0;
    Real hpwl = 0.0;
    Real overflow = 0.0;
};

struct RecoveryStats {
    int sweeps = 0;
    int moves = 0;
    int objective_evaluations = 0;
    int net_block_attempts = 0;
    int net_block_moves = 0;
    int net_block_nodes = 0;
    int net_block_contraction_attempts = 0;
    int net_block_contraction_moves = 0;
    int net_block_contraction_nodes = 0;
    int compact_direction_attempts = 0;
    int compact_direction_moves = 0;
    Real initial_hpwl = 0.0;
    Real initial_overflow = 0.0;
    Real final_hpwl = 0.0;
    Real final_overflow = 0.0;
    double wall_seconds = 0.0;
    std::vector<RecoverySweepStats> trajectory;
};

RecoveryStats recover_hpwl_under_overflow(
    Database& db, ExactOverlapDensity& density, const RecoveryConfig& config);

}  // namespace ea
