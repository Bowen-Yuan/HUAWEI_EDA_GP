#pragma once

#include "types.h"

#include <vector>

namespace dpcpp {

Real exact_hpwl(const Database& db, int degree_limit = -1);
Real exact_hpwl_subgradient(const Database& db, int gradient_degree_limit,
                            std::vector<Real>* grad_x,
                            std::vector<Real>* grad_y);
Real weighted_average_wirelength(const Database& db, Real gamma,
                                 int degree_limit,
                                 std::vector<Real>* grad_x,
                                 std::vector<Real>* grad_y);

}  // namespace dpcpp
