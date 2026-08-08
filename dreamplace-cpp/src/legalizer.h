#pragma once

#include "types.h"

#include <string>
#include <vector>

namespace dpcpp {

struct LegalizeConfig {
    int row_search_limit = 64;
    int detailed_passes = 3;
    int k_reorder_size = 4;
    int global_swap_passes = 1;
    int independent_set_size = 4;
    int detailed_outer_rounds = 1;
    int cell_insertion_passes = 0;
    int cell_insertion_window = 16;
    int projected_subgradient_passes = 0;
    Real projected_step_sites = 4.0;
    int constrained_bundle_passes = 0;
    int constrained_bundle_size = 8;
    Real constrained_bundle_step_sites = 4.0;
    int row_relegalization_passes = 0;
    bool use_hungarian_matching = false;
    Real minimum_relative_improvement = 1.0e-6;
    bool run_abacus = true;
    bool run_detailed = true;
    bool run_dreamplace_detailed = false;
    std::string snapshot_dir;
};

struct LegalizationProxy {
    Real mean_row_distance = 0.0;
    Real segment_overflow = 0.0;
};

struct LegalityResult {
    bool legal = false;
    int boundary_errors = 0;
    int alignment_errors = 0;
    int overlap_errors = 0;
    std::string first_error;
};

struct LegalizeResult {
    Real hpwl_before = 0.0;
    Real hpwl_after_greedy = 0.0;
    Real hpwl_after_abacus = 0.0;
    Real hpwl_after_k_reorder = 0.0;
    Real hpwl_after_global_swap = 0.0;
    Real hpwl_after_independent_set = 0.0;
    Real hpwl_after_constrained_bundle = 0.0;
    Real hpwl_after_detailed = 0.0;
    LegalityResult legality;
};

LegalityResult check_legality(const Database& db, Real tolerance = 1.0e-5);
LegalizationProxy evaluate_legalization_proxy(const Database& db);
void compute_legalization_force(const Database& db, Real congestion_gain,
                                std::vector<Real>& grad_x,
                                std::vector<Real>& grad_y);
LegalizeResult legalize_and_refine(Database& db, const LegalizeConfig& config);

}  // namespace dpcpp
