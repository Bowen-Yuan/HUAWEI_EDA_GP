#include "epsilon_active/step_policy.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ea {
namespace {

class ConstantController final : public StepController {
public:
    explicit ConstantController(const StepControllerConfig& config)
        : config_(config) {}

    void reset(int, Real bin_size) override {
        maximum_delta_ = config_.maximum_delta_bins * bin_size;
    }

    StepDecision propose(int) override {
        return {config_.learning_rate, maximum_delta_};
    }

    void observe(const StepObservation&) override {}

private:
    StepControllerConfig config_;
    Real maximum_delta_ = 0.0;
};

class CosineDecayController final : public StepController {
public:
    explicit CosineDecayController(const StepControllerConfig& config)
        : config_(config) {}

    void reset(int total_iterations, Real bin_size) override {
        if (total_iterations < 1) {
            throw std::invalid_argument("cosine step policy needs a positive iteration budget");
        }
        total_iterations_ = total_iterations;
        maximum_delta_ = config_.maximum_delta_bins * bin_size;
    }

    StepDecision propose(int iteration) override {
        const Real q = std::clamp(
            static_cast<Real>(iteration) / static_cast<Real>(total_iterations_),
            static_cast<Real>(0.0), static_cast<Real>(1.0));
        const Real lr_min = config_.learning_rate_min_ratio * config_.learning_rate;
        const Real lr = lr_min + 0.5 * (config_.learning_rate - lr_min) *
                                    (1.0 + std::cos(static_cast<Real>(M_PI) * q));
        return {lr, maximum_delta_};
    }

    void observe(const StepObservation&) override {}

private:
    StepControllerConfig config_;
    int total_iterations_ = 0;
    Real maximum_delta_ = 0.0;
};

class TrustRadiusController final : public StepController {
public:
    explicit TrustRadiusController(const StepControllerConfig& config)
        : config_(config) {}

    void reset(int, Real bin_size) override {
        bin_size_ = bin_size;
        radius_ = config_.trust_radius_bins * bin_size;
        clamp_radius();
        accept_streak_ = 0;
    }

    StepDecision propose(int) override { return {config_.learning_rate, radius_}; }

    void observe(const StepObservation& observation) override {
        if (observation.accepted) {
            ++accept_streak_;
            if (accept_streak_ >= config_.grow_after_accepts) {
                radius_ *= config_.grow_factor;
                accept_streak_ = 0;
                clamp_radius();
            }
        } else {
            radius_ *= config_.shrink_factor;
            accept_streak_ = 0;
            clamp_radius();
        }
    }

    void reset_streaks() override { accept_streak_ = 0; }

private:
    void clamp_radius() {
        radius_ = std::clamp(radius_,
                             config_.radius_min_bins * bin_size_,
                             config_.radius_max_bins * bin_size_);
    }

    StepControllerConfig config_;
    Real bin_size_ = 1.0;
    Real radius_ = 0.0;
    int accept_streak_ = 0;
};

}  // namespace

StepPolicyKind parse_step_policy(const std::string& name) {
    if (name == "constant") return StepPolicyKind::Constant;
    if (name == "cosine") return StepPolicyKind::CosineDecay;
    if (name == "trust") return StepPolicyKind::TrustRadius;
    throw std::invalid_argument("unknown step policy: " + name);
}

const char* step_policy_name(StepPolicyKind kind) noexcept {
    switch (kind) {
    case StepPolicyKind::Constant: return "constant";
    case StepPolicyKind::CosineDecay: return "cosine";
    case StepPolicyKind::TrustRadius: return "trust";
    }
    return "unknown";
}

void StepControllerConfig::validate() const {
    if (!(learning_rate > 0.0) || !(maximum_delta_bins > 0.0)) {
        throw std::invalid_argument("invalid step policy scale parameters");
    }
    if (!(learning_rate_min_ratio > 0.0) || !(learning_rate_min_ratio < 1.0)) {
        throw std::invalid_argument("learning_rate_min_ratio must be in (0, 1)");
    }
    if (!(trust_radius_bins > 0.0) || !(radius_min_bins > 0.0) ||
        !(radius_max_bins >= radius_min_bins)) {
        throw std::invalid_argument("invalid trust radius bounds");
    }
    if (!(grow_factor > 1.0) || !(shrink_factor > 0.0) || !(shrink_factor < 1.0)) {
        throw std::invalid_argument("invalid trust radius growth factors");
    }
    if (grow_after_accepts < 1) {
        throw std::invalid_argument("grow_after_accepts must be positive");
    }
}

std::unique_ptr<StepController> make_step_controller(
    StepPolicyKind kind, const StepControllerConfig& config) {
    config.validate();
    switch (kind) {
    case StepPolicyKind::Constant:
        return std::make_unique<ConstantController>(config);
    case StepPolicyKind::CosineDecay:
        return std::make_unique<CosineDecayController>(config);
    case StepPolicyKind::TrustRadius:
        return std::make_unique<TrustRadiusController>(config);
    }
    throw std::invalid_argument("unknown step policy kind");
}

}  // namespace ea
