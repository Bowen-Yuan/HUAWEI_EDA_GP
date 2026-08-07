#pragma once

#include "types.h"

#include <string>

namespace dpcpp {

struct LegalizeConfig {
    int row_search_limit = 64;
    int detailed_passes = 3;
    int k_reorder_size = 4;
    int global_swap_passes = 1;
    int independent_set_size = 4;
    bool run_abacus = true;
    bool run_detailed = true;
    bool run_dreamplace_detailed = false;
    std::string snapshot_dir;
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
    Real hpwl_after_detailed = 0.0;
    LegalityResult legality;
};

LegalityResult check_legality(const Database& db, Real tolerance = 1.0e-5);
LegalizeResult legalize_and_refine(Database& db, const LegalizeConfig& config);

}  // namespace dpcpp
