#include "epsilon_active/lambda_controller.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ea {

LambdaPolicy parse_lambda_policy(const std::string& name) {
    if (name == "dreamplace") return LambdaPolicy::Dreamplace;
    if (name == "trajectory") return LambdaPolicy::Trajectory;
    if (name == "ratio") return LambdaPolicy::Ratio;
    throw std::invalid_argument("unknown lambda policy: " + name);
}

const char* lambda_policy_name(LambdaPolicy policy) noexcept {
    switch (policy) {
    case LambdaPolicy::Dreamplace: return "dreamplace";
    case LambdaPolicy::Trajectory: return "trajectory";
    case LambdaPolicy::Ratio: return "ratio";
    }
    return "unknown";
}

LambdaController::LambdaController(LambdaConfig config) : config_(config) {}

void LambdaController::initialize(const Database& db,
                                  const std::vector<Real>& wire_x,
                                  const std::vector<Real>& wire_y,
                                  const std::vector<Real>& density_x,
                                  const std::vector<Real>& density_y,
                                  Real hpwl, Real overflow) {
    Real wire_l1 = 0.0;
    Real density_l1 = 0.0;
    #pragma omp parallel for reduction(+:wire_l1,density_l1) schedule(static)
    for (int movable = 0; movable < static_cast<int>(db.movable_ids.size()); ++movable) {
        const int id = db.movable_ids[movable];
        wire_l1 += std::abs(wire_x[id]) + std::abs(wire_y[id]);
        density_l1 += std::abs(density_x[id]) + std::abs(density_y[id]);
    }
    base_ = config_.density_weight_scale * wire_l1 /
            std::max<Real>(density_l1, 1.0e-30);
    if (config_.initial_effective > 0.0 && base_ > 0.0) {
        control_ = std::clamp(config_.initial_effective / base_,
                              config_.control_min, config_.control_max);
    } else {
        control_ = std::clamp(config_.initial_control,
                              config_.control_min, config_.control_max);
    }
    previous_hpwl_ = hpwl;
    previous_overflow_ = overflow;
    initial_overflow_ = overflow;
    integral_error_ = 0.0;
}

void LambdaController::update(int density_step, Real hpwl, Real overflow) {
    if (config_.policy == LambdaPolicy::Dreamplace) {
        const Real delta = hpwl - previous_hpwl_;
        Real multiplier = 1.0;
        if (delta < 0.0) {
            multiplier = 1.05 * std::max(std::pow(0.9999, density_step), 0.98);
        } else {
            multiplier = std::clamp(
                std::pow(1.05, 1.0 - delta / config_.reference_hpwl_delta),
                0.95, 1.05);
        }
        control_ = std::clamp(control_ * multiplier,
                              config_.control_min, config_.control_max);
    } else if (config_.update_interval > 0 &&
               density_step % config_.update_interval == 0) {
        if (config_.policy == LambdaPolicy::Ratio) {
            const Real hpwl_target = config_.hpwl_baseline > 0.0
                ? config_.hpwl_baseline : hpwl;
            const Real hpwl_ratio = hpwl / std::max<Real>(hpwl_target, 1.0);
            const Real overflow_ratio = overflow /
                std::max<Real>(config_.overflow_baseline, 1.0e-9);
            Real factor = 1.0;
            if (overflow_ratio > 1.05 * hpwl_ratio) factor = 1.5;
            else if (hpwl_ratio > 1.05 * overflow_ratio) factor = 2.0 / 3.0;
            if (overflow > 1.05 * config_.stop_overflow) factor = std::max(factor, 1.25);
            control_ = std::clamp(control_ * factor,
                                  config_.control_min, config_.control_max);
        } else {
            const Real progress = std::min<Real>(
                1.0, density_step /
                static_cast<Real>(std::max(1, config_.trajectory_horizon)));
            const Real schedule = progress * progress * (3.0 - 2.0 * progress);
            const Real desired = config_.stop_overflow +
                (initial_overflow_ - config_.stop_overflow) * (1.0 - schedule);
            const Real scale = std::max<Real>(config_.stop_overflow, 1.0e-4);
            const Real error = std::clamp((overflow - desired) / scale, -2.0, 4.0);
            const Real desired_drop = (previous_overflow_ - desired) / scale;
            const Real actual_drop = (previous_overflow_ - overflow) / scale;
            const Real trend_error = std::clamp(desired_drop - actual_drop, -2.0, 2.0);
            if (progress < 1.0) {
                integral_error_ = std::clamp(integral_error_ + error, -3.0, 3.0);
            } else {
                integral_error_ = 0.0;
            }
            Real per_iteration = 0.045 + 0.08 * error +
                0.003 * integral_error_ + 0.04 * trend_error;
            if (overflow <= config_.stop_overflow && config_.hpwl_baseline > 0.0) {
                per_iteration -= 0.01 * std::max<Real>(
                    0.0, hpwl / config_.hpwl_baseline - 1.0);
            }
            per_iteration = std::clamp(per_iteration, -0.05, 0.06);
            const Real log_delta = std::clamp(
                per_iteration * config_.update_interval, -1.0, 1.25);
            control_ = std::clamp(control_ * std::exp(log_delta),
                                  config_.control_min, config_.control_max);
        }
        previous_overflow_ = overflow;
    }
    previous_hpwl_ = hpwl;
}

}  // namespace ea
