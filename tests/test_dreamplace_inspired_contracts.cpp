#include "epsilon_active/spectral_step.hpp"
#include "epsilon_active/epsilon_continuation.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool nearly_equal(ea::Real a, ea::Real b, ea::Real tol) {
    return std::abs(a - b) <= tol;
}

ea::SpectralStepDecision bb_update(ea::SpectralStep& step,
                                   const std::vector<ea::Real>& s,
                                   const std::vector<ea::Real>& y) {
    return step.update(s, y);
}

}  // namespace

int main() {
    // BB1: positive curvature formula alpha_bb1 = sTs / sTy.
    {
        ea::SpectralStepConfig c;
        c.method = "bb1";
        c.base_learning_rate = 1.0;
        c.alpha_min_ratio = 0.125;
        c.alpha_max_ratio = 8.0;
        c.max_step_growth = 1.0e9;
        ea::SpectralStep step(c);
        // s = 2u, y = u with unit u: sTs=4, sTy=2, yTy=2 -> bb1 = 2.
        std::vector<ea::Real> s{2.0, 0.0, 2.0, 0.0};
        std::vector<ea::Real> y{1.0, 0.0, 1.0, 0.0};
        const auto d = bb_update(step, s, y);
        if (!d.valid) return 1;
        if (!nearly_equal(d.alpha_bb1, 2.0, 1e-12)) return 2;
        if (!nearly_equal(d.alpha_used, 2.0, 1e-12)) return 3;
    }
    // BB2: negative curvature falls back to base and flags a restart.
    {
        ea::SpectralStepConfig c;
        c.method = "bb2";
        c.base_learning_rate = 0.5;
        ea::SpectralStep step(c);
        std::vector<ea::Real> s{1.0, 0.0};
        std::vector<ea::Real> y{-1.0, 0.0};  // sTy = -1 < 0
        const auto d = bb_update(step, s, y);
        if (d.valid) return 4;
        if (!d.restarted) return 5;
        if (!nearly_equal(d.alpha_used, 0.5, 1e-12)) return 6;
    }
    // Zero y must never produce NaN/Inf.
    {
        ea::SpectralStepConfig c;
        ea::SpectralStep step(c);
        std::vector<ea::Real> s{1.0, 1.0};
        std::vector<ea::Real> y(2, 0.0);
        const auto d = bb_update(step, s, y);
        if (d.valid) return 7;
        if (!std::isfinite(d.alpha_used) || !std::isfinite(d.sTy) ||
            !std::isfinite(d.yTy)) return 8;
    }
    // Clipping keeps alpha inside [min,max] bounds.
    {
        ea::SpectralStepConfig c;
        c.method = "bb2";
        c.base_learning_rate = 1.0;
        c.alpha_min_ratio = 0.125;
        c.alpha_max_ratio = 2.0;
        c.max_step_growth = 1.0e9;
        ea::SpectralStep step(c);
        // bb2 = sTy/yTy = 100/1 = 100 -> clipped to 2.
        std::vector<ea::Real> s{10.0, 0.0};
        std::vector<ea::Real> y{0.1, 0.0};  // sTy = 1, yTy = 0.01
        const auto d = bb_update(step, s, y);
        if (!nearly_equal(d.alpha_used, 2.0, 1e-12)) return 9;
        // tiny alpha clipped to 0.125
        std::vector<ea::Real> s2{0.1, 0.0};
        std::vector<ea::Real> y2{10.0, 0.0};
        const auto d2 = bb_update(step, s2, y2);
        if (!nearly_equal(d2.alpha_used, 0.125, 1e-12)) return 10;
    }
    // Growth limiter: consecutive alpha cannot exceed max_step_growth.
    {
        ea::SpectralStepConfig c;
        c.method = "bb2";
        c.base_learning_rate = 1.0;
        c.max_step_growth = 2.0;
        ea::SpectralStep step(c);
        std::vector<ea::Real> s{1.0, 0.0};
        std::vector<ea::Real> y{1.0, 0.0};  // bb2 = 1
        const auto d1 = bb_update(step, s, y);
        if (!nearly_equal(d1.alpha_used, 1.0, 1e-12)) return 11;
        std::vector<ea::Real> s2{10.0, 0.0};
        std::vector<ea::Real> y2{1.0, 0.0};  // bb2 = 10, capped to 2
        const auto d2 = bb_update(step, s2, y2);
        if (!nearly_equal(d2.alpha_used, 2.0, 1e-12)) return 12;
    }
    // BB7: the estimator is a pure scalar computation; running it must not
    // mutate its inputs (no rollback / no layout semantics involved).
    {
        ea::SpectralStepConfig c;
        ea::SpectralStep step(c);
        std::vector<ea::Real> s{1.0, 2.0, 3.0};
        std::vector<ea::Real> y{1.0, 1.0, 1.0};
        const auto d = bb_update(step, s, y);
        (void)d;
        if (s[0] != 1.0 || s[2] != 3.0 || y[1] != 1.0) return 13;
    }

    // EC-C1/C2: monotonicity and endpoints of the epsilon schedule.
    {
        ea::EpsilonContinuationConfig c;
        c.overflow_low = 0.07;
        c.overflow_high = 0.15;
        c.wire_min_bins = 0.25;
        c.wire_max_bins = 1.0;
        ea::EpsilonContinuation ec(c);
        const auto low = ec.observe(0.05);
        const auto mid = ec.observe(0.11);
        const auto high = ec.observe(0.20);
        if (!nearly_equal(low.wire_epsilon_bins, 0.25, 1e-12)) return 14;
        if (!nearly_equal(high.wire_epsilon_bins, 1.0, 1e-12)) return 15;
        if (!(low.wire_epsilon_bins <= mid.wire_epsilon_bins &&
              mid.wire_epsilon_bins <= high.wire_epsilon_bins)) return 16;
        if (!(low.density_epsilon_bins <= mid.density_epsilon_bins &&
              mid.density_epsilon_bins <= high.density_epsilon_bins)) return 17;
        // q mid with smoothstep: q = (0.11-0.07)/0.08 = 0.5 -> h = 0.5
        if (!nearly_equal(mid.wire_epsilon_bins, 0.625, 1e-12)) return 18;
    }
    // EC-C4: band change triggers restart flag only when crossing bands.
    {
        ea::EpsilonContinuationConfig c;
        c.restart_on_band_change = true;
        ea::EpsilonContinuation ec(c);
        const auto a = ec.observe(0.10);  // q=0.375 -> medium band
        if (a.restart) return 19;         // first observation: no restart
        const auto b = ec.observe(0.11);  // same band
        if (b.restart) return 20;
        const auto d = ec.observe(0.14);  // q=0.875 -> wide band
        if (!d.restart) return 21;
    }

    std::cout << "dreamplace-inspired contracts passed\n";
}

