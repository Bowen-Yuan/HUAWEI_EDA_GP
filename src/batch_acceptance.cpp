#include "epsilon_active/batch_acceptance.hpp"

#include <stdexcept>

namespace ea {

BatchAcceptanceResult accept_exact_batch(
    const std::vector<Real>& positions,
    const std::vector<Real>& delta,
    Real current_hpwl,
    Real current_overflow,
    const BatchAcceptanceConfig& config,
    const PositionApplier& apply,
    const ExactBatchEvaluator& evaluate) {
    if (positions.size() != delta.size()) {
        throw std::invalid_argument("batch delta has wrong size");
    }
    if (config.max_trials <= 0 || config.shrink <= 0.0 ||
        config.shrink >= 1.0 || config.overflow_cap < 0.0) {
        throw std::invalid_argument("invalid batch acceptance configuration");
    }

    const std::vector<Real> baseline = positions;
    BatchAcceptanceResult result;
    result.hpwl = current_hpwl;
    result.overflow = current_overflow;
    std::vector<Real> candidate(positions.size());
    Real scale = 1.0;
    for (int trial = 0; trial < config.max_trials; ++trial) {
        for (std::size_t coordinate = 0; coordinate < positions.size(); ++coordinate) {
            candidate[coordinate] = baseline[coordinate] - scale * delta[coordinate];
        }
        apply(candidate);
        const auto [hpwl, overflow] = evaluate();
        ++result.trials;
        const bool feasible_descent =
            overflow <= config.overflow_cap && hpwl < current_hpwl;
        const bool infeasible_progress =
            config.allow_infeasible && current_overflow > config.overflow_cap &&
            overflow < current_overflow &&
            hpwl <= current_hpwl * (1.0 + config.max_infeasible_hpwl_increase);
        if (feasible_descent || infeasible_progress) {
            result.accepted = true;
            result.scale = scale;
            result.hpwl = hpwl;
            result.overflow = overflow;
            return result;
        }
        scale *= config.shrink;
    }

    apply(baseline);
    return result;
}

}  // namespace ea
