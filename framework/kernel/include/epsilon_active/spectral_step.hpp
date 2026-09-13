#pragma once

#include "epsilon_active/types.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ea {

// Spectral (Barzilai-Borwein) step-size estimation for non-smooth global
// placement.  This is a pure step-size adaptation: it changes only the
// scalar learning rate fed to the optimizer, never the direction oracle,
// the exact objectives, or the acceptance contract.  Invalid curvature
// falls back to the base learning rate (a state restart, not a rollback).
struct SpectralStepConfig {
    std::string method = "bb2";  // "bb1" or "bb2"
    Real base_learning_rate = 0.002;
    Real alpha_min_ratio = 0.125;
    Real alpha_max_ratio = 8.0;
    Real max_step_growth = 2.0;
    Real curvature_epsilon = 1.0e-14;
    Real lipschitz_factor = 1.0;
    Real raw_sanity_max = 1.0e30;
};

struct SpectralStepDecision {
    bool valid = false;
    bool restarted = false;
    Real alpha_bb1 = 0.0;
    Real alpha_bb2 = 0.0;
    Real alpha_lip = 0.0;
    Real alpha_used = 0.0;
    Real sTy = 0.0;
    Real yTy = 0.0;
};

class SpectralStep {
public:
    explicit SpectralStep(SpectralStepConfig config)
        : config_(std::move(config)),
          previous_alpha_(config_.base_learning_rate) {
        if (config_.method != "bb1" && config_.method != "bb2")
            throw std::invalid_argument("spectral step requires method bb1 or bb2");
        if (config_.base_learning_rate <= 0.0 || config_.alpha_min_ratio <= 0.0 ||
            config_.alpha_max_ratio < config_.alpha_min_ratio ||
            config_.max_step_growth < 1.0 || config_.curvature_epsilon < 0.0)
            throw std::invalid_argument("invalid spectral step configuration");
    }

    // s = x_k - x_{k-1}, y = g_k - g_{k-1} (both flat [x..., y...] vectors).
    SpectralStepDecision update(const std::vector<Real>& s,
                                const std::vector<Real>& y) {
        SpectralStepDecision decision;
        decision.alpha_used = config_.base_learning_rate;
        if (s.size() != y.size() || s.empty()) {
            decision.restarted = restart_invalid_state();
            return decision;
        }
        Real sTy = 0.0, sTs = 0.0, yTy = 0.0, s_norm = 0.0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            sTy += s[i] * y[i];
            sTs += s[i] * s[i];
            yTy += y[i] * y[i];
        }
        s_norm = std::sqrt(sTs);
        const Real y_norm = std::sqrt(yTy);
        decision.sTy = sTy;
        decision.yTy = yTy;
        decision.alpha_bb1 = sTs / std::max<Real>(sTy, 1.0e-300);
        decision.alpha_bb2 = sTy / std::max<Real>(yTy, 1.0e-300);
        decision.alpha_lip = s_norm / std::max<Real>(y_norm, 1.0e-300);

        const bool curvature_ok =
            std::isfinite(sTy) && std::isfinite(yTy) &&
            sTy > config_.curvature_epsilon && yTy > config_.curvature_epsilon;
        if (!curvature_ok) {
            decision.restarted = restart_invalid_state();
            return decision;
        }

        Real raw = config_.method == "bb1" ? decision.alpha_bb1
                                           : decision.alpha_bb2;
        raw = std::min(raw, config_.lipschitz_factor * decision.alpha_lip);
        const bool sane = std::isfinite(raw) && raw > 0.0 &&
                          raw <= config_.raw_sanity_max;
        if (!sane) {
            decision.restarted = restart_invalid_state();
            return decision;
        }

        const Real lo = config_.alpha_min_ratio * config_.base_learning_rate;
        const Real hi = config_.alpha_max_ratio * config_.base_learning_rate;
        Real alpha = std::clamp(raw, lo, hi);
        // Growth limiter against sudden active-set switches.
        const Real growth_cap = config_.max_step_growth * previous_alpha_;
        alpha = std::min(alpha, growth_cap);
        alpha = std::clamp(alpha, lo, hi);
        previous_alpha_ = alpha;
        decision.valid = true;
        decision.alpha_used = alpha;
        return decision;
    }

    Real previous_alpha() const noexcept { return previous_alpha_; }

private:
    bool restart_invalid_state() {
        previous_alpha_ = config_.base_learning_rate;
        return true;
    }

    SpectralStepConfig config_;
    // Seeded with the base rate so the growth limiter never collapses the
    // first accepted step to zero.
    Real previous_alpha_;
};

}  // namespace ea
