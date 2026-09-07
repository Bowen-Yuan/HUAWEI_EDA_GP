#pragma once

#include "epsilon_active/types.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ea {

// Keep all global-view update rules behind this one interface.  The lab
// runner deliberately does not grow a second optimizer hierarchy.
enum class OptimizerKind {
    Adam, AMSGrad, AdaGrad, HeavyBall, SGD, NormalizedSGD, DualAveraging
};

OptimizerKind parse_optimizer(const std::string& name);
const char* optimizer_name(OptimizerKind kind) noexcept;

class Optimizer {
public:
    virtual ~Optimizer() = default;
    virtual const char* name() const noexcept = 0;
    virtual void reset(std::size_t dimensions) = 0;
    virtual void compute_delta(const std::vector<Real>& gradient,
                               Real learning_rate, Real maximum_delta,
                               std::vector<Real>& delta) = 0;
};

std::unique_ptr<Optimizer> make_optimizer(OptimizerKind kind, Real beta1,
                                          Real beta2, Real momentum,
                                          Real numerical_epsilon);

}  // namespace ea
