#include "epsilon_active/optimizer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ea {
namespace {

class StatefulOptimizer final : public Optimizer {
public:
    StatefulOptimizer(OptimizerKind kind, Real beta1, Real beta2,
                      Real momentum, Real numerical_epsilon)
        : kind_(kind), beta1_(beta1), beta2_(beta2), momentum_(momentum),
          numerical_epsilon_(numerical_epsilon) {}

    const char* name() const noexcept override { return optimizer_name(kind_); }

    void reset(std::size_t dimensions) override {
        first_.assign(dimensions, 0.0);
        second_.assign(dimensions, 0.0);
        maximum_second_.assign(dimensions, 0.0);
        age_ = 0;
    }

    void compute_delta(const std::vector<Real>& gradient, Real learning_rate,
                       Real maximum_delta, std::vector<Real>& delta) override {
        if (first_.size() != gradient.size()) reset(gradient.size());
        delta.resize(gradient.size());
        ++age_;
        const Real first_bias = std::max<Real>(1.0e-12, 1.0 - std::pow(beta1_, age_));
        const Real second_bias = std::max<Real>(1.0e-12, 1.0 - std::pow(beta2_, age_));
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(gradient.size()); ++i) {
            const Real g = gradient[i];
            Real value = 0.0;
            switch (kind_) {
            case OptimizerKind::SGD:
                value = learning_rate * g;
                break;
            case OptimizerKind::HeavyBall:
                first_[i] = momentum_ * first_[i] + (1.0 - momentum_) * g;
                value = learning_rate * first_[i];
                break;
            case OptimizerKind::AdaGrad:
                second_[i] += g * g;
                value = learning_rate * g /
                        (std::sqrt(second_[i]) + numerical_epsilon_);
                break;
            case OptimizerKind::Adam:
            case OptimizerKind::AMSGrad: {
                first_[i] = beta1_ * first_[i] + (1.0 - beta1_) * g;
                second_[i] = beta2_ * second_[i] + (1.0 - beta2_) * g * g;
                const Real first_hat = first_[i] / first_bias;
                const Real second_hat = second_[i] / second_bias;
                Real denominator = second_hat;
                if (kind_ == OptimizerKind::AMSGrad) {
                    maximum_second_[i] = std::max(maximum_second_[i], second_hat);
                    denominator = maximum_second_[i];
                }
                value = learning_rate * first_hat /
                        (std::sqrt(denominator) + numerical_epsilon_);
                break;
            }
            }
            delta[i] = std::clamp(value, -maximum_delta, maximum_delta);
        }
    }

private:
    OptimizerKind kind_;
    Real beta1_;
    Real beta2_;
    Real momentum_;
    Real numerical_epsilon_;
    int age_ = 0;
    std::vector<Real> first_;
    std::vector<Real> second_;
    std::vector<Real> maximum_second_;
};

}  // namespace

OptimizerKind parse_optimizer(const std::string& name) {
    if (name == "adam") return OptimizerKind::Adam;
    if (name == "amsgrad") return OptimizerKind::AMSGrad;
    if (name == "adagrad") return OptimizerKind::AdaGrad;
    if (name == "heavy-ball") return OptimizerKind::HeavyBall;
    if (name == "sgd") return OptimizerKind::SGD;
    throw std::invalid_argument("unknown optimizer: " + name);
}

const char* optimizer_name(OptimizerKind kind) noexcept {
    switch (kind) {
    case OptimizerKind::Adam: return "adam";
    case OptimizerKind::AMSGrad: return "amsgrad";
    case OptimizerKind::AdaGrad: return "adagrad";
    case OptimizerKind::HeavyBall: return "heavy-ball";
    case OptimizerKind::SGD: return "sgd";
    }
    return "unknown";
}

std::unique_ptr<Optimizer> make_optimizer(OptimizerKind kind, Real beta1,
                                          Real beta2, Real momentum,
                                          Real numerical_epsilon) {
    if (beta1 < 0.0 || beta1 >= 1.0 || beta2 < 0.0 || beta2 >= 1.0 ||
        momentum < 0.0 || momentum >= 1.0 || numerical_epsilon <= 0.0) {
        throw std::invalid_argument("invalid optimizer hyperparameters");
    }
    return std::make_unique<StatefulOptimizer>(
        kind, beta1, beta2, momentum, numerical_epsilon);
}

}  // namespace ea

