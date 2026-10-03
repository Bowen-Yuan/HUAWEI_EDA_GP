#pragma once

#include "epsilon_active/types.hpp"

#include <functional>
#include <utility>
#include <vector>

namespace ea {

struct BatchAcceptanceConfig {
    bool enabled = false;
    // Allows exact line-search progress before the 7% cap is reached. The
    // accepted candidate must reduce exact overflow and stay within the
    // configured relative HPWL damage budget.
    bool allow_infeasible = false;
    Real max_infeasible_hpwl_increase = 0.01;
    int max_trials = 12;
    Real shrink = 0.5;
    Real overflow_cap = 0.07;
};

struct BatchAcceptanceResult {
    bool accepted = false;
    Real scale = 0.0;
    int trials = 0;
    Real hpwl = 0.0;
    Real overflow = 0.0;
};

using PositionApplier = std::function<void(const std::vector<Real>&)>;
using ExactBatchEvaluator = std::function<std::pair<Real, Real>()>;

BatchAcceptanceResult accept_exact_batch(
    const std::vector<Real>& positions,
    const std::vector<Real>& delta,
    Real current_hpwl,
    Real current_overflow,
    const BatchAcceptanceConfig& config,
    const PositionApplier& apply,
    const ExactBatchEvaluator& evaluate);

}  // namespace ea
