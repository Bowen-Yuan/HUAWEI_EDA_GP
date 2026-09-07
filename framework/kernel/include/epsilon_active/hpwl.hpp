#pragma once

#include "epsilon_active/types.hpp"

#include <vector>

namespace ea {

class ExactHpwl {
public:
    explicit ExactHpwl(const Database& db);

    Real evaluate(Real epsilon, Real active_power, int degree_limit,
                  std::vector<Real>* grad_x, std::vector<Real>* grad_y);

private:
    const Database& db_;
    std::vector<Real> pin_grad_x_;
    std::vector<Real> pin_grad_y_;
};

}  // namespace ea

