#pragma once

#include "epsilon_active/types.hpp"

#include <vector>

namespace ea {

// Direction handling for nets whose degree exceeds the configured limit.
// Ignore keeps the historical behavior (no direction contribution).
// ExactExtrema contributes the epsilon=0 exact extremal-face subgradient of
// those nets: min-side ties split -w_e, max-side ties split +w_e.  This is a
// legitimate exact non-smooth subgradient and never a smooth surrogate.
enum class HighDegreeMode {
    Ignore,
    ExactExtrema
};

class ExactHpwl {
public:
    explicit ExactHpwl(const Database& db);

    Real evaluate(Real epsilon, Real active_power, int degree_limit,
                  std::vector<Real>* grad_x, std::vector<Real>* grad_y);
    Real evaluate(Real epsilon, Real active_power, int degree_limit,
                  HighDegreeMode high_degree_mode,
                  std::vector<Real>* grad_x, std::vector<Real>* grad_y);

private:
    const Database& db_;
    std::vector<Real> pin_grad_x_;
    std::vector<Real> pin_grad_y_;
};

}  // namespace ea
