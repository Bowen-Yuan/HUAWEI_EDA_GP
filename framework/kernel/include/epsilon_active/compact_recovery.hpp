#pragma once

#include "epsilon_active/density.hpp"
#include "epsilon_active/types.hpp"

#include <filesystem>
#include <vector>

namespace ea {

struct CompactRecoveryConfig {
    int sweeps = 0;
    int line_search_steps = 6;
    Real max_step_bins = 32.0;
    Real overflow_cap = -1.0;
    Real hpwl_tolerance = 1.0e-10;
    bool axis_separated = true;
    bool fixed_terminal_aware = true;
    int snapshot_every = 0;
    std::filesystem::path output_dir;
};

struct CompactRecoverySweepStats {
    int sweep = 0;
    int moves = 0;
    Real hpwl = 0.0;
    Real overflow = 0.0;
    Real compactness = 0.0;
};

struct CompactRecoveryStats {
    int sweeps = 0;
    int moves = 0;
    int candidates = 0;
    int accepted_candidates = 0;
    int objective_evaluations = 0;
    Real initial_hpwl = 0.0;
    Real initial_overflow = 0.0;
    Real final_hpwl = 0.0;
    Real final_overflow = 0.0;
    Real initial_compactness = 0.0;
    Real final_compactness = 0.0;
    double wall_seconds = 0.0;
    std::vector<CompactRecoverySweepStats> trajectory;
};

CompactRecoveryStats compact_support_contraction(
    Database& db, ExactOverlapDensity& density,
    const CompactRecoveryConfig& config);

}  // namespace ea
