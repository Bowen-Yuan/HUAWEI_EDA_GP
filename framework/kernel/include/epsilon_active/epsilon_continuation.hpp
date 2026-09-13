#pragma once

#include "epsilon_active/types.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ea {

// Overflow-driven epsilon-active radius continuation (DREAMPlace-inspired).
// Only the direction oracle radius changes; exact HPWL/overflow metrics and
// the final audit always use epsilon = 0 and are unaffected.  High overflow
// widens the active face; approaching the target tightens it back.
struct EpsilonContinuationConfig {
    Real overflow_low = 0.07;    // ratio, not percent
    Real overflow_high = 0.15;   // ratio, not percent
    Real wire_min_bins = 0.25;
    Real wire_max_bins = 1.0;
    Real density_min_bins = 0.0;
    Real density_max_bins = 0.5;
    bool restart_on_band_change = true;
};

struct EpsilonDecision {
    Real q = 0.0;                     // normalized schedule coordinate [0,1]
    Real wire_epsilon_bins = 0.0;
    Real density_epsilon_bins = 0.0;
    int band = 0;                     // 0 narrow, 1 medium, 2 wide
    bool restart = false;             // band changed this observation
};

class EpsilonContinuation {
public:
    explicit EpsilonContinuation(EpsilonContinuationConfig config)
        : config_(config) {
        if (!(config_.overflow_low < config_.overflow_high) ||
            config_.wire_min_bins < 0.0 || config_.density_min_bins < 0.0 ||
            config_.wire_max_bins < config_.wire_min_bins ||
            config_.density_max_bins < config_.density_min_bins)
            throw std::invalid_argument("invalid epsilon continuation configuration");
    }

    static Real smoothstep(Real q) {
        return q * q * (3.0 - 2.0 * q);
    }

    // Observes the current exact overflow ratio and produces the epsilon
    // radii (in bin units) plus the band-change restart flag.
    EpsilonDecision observe(Real overflow_ratio) {
        const Real span = std::max<Real>(
            config_.overflow_high - config_.overflow_low, 1.0e-12);
        const Real q = std::clamp(
            (overflow_ratio - config_.overflow_low) / span, 0.0, 1.0);
        EpsilonDecision decision;
        decision.q = q;
        const Real h = smoothstep(q);
        decision.wire_epsilon_bins = config_.wire_min_bins +
            h * (config_.wire_max_bins - config_.wire_min_bins);
        decision.density_epsilon_bins = config_.density_min_bins +
            h * (config_.density_max_bins - config_.density_min_bins);
        decision.band = q > 0.66 ? 2 : (q > 0.33 ? 1 : 0);
        if (has_previous_ && decision.band != previous_band_) {
            decision.restart = config_.restart_on_band_change;
        }
        previous_band_ = decision.band;
        has_previous_ = true;
        return decision;
    }

private:
    EpsilonContinuationConfig config_;
    int previous_band_ = 0;
    bool has_previous_ = false;
};

}  // namespace ea

