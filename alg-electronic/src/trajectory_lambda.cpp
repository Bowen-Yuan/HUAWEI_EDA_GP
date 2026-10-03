#include "trajectory_lambda.h"

#include <algorithm>
#include <cmath>

TrajectoryLambdaController::TrajectoryLambdaController(
    TrajectoryLambdaConfig config) : config_(config) {}

TrajectoryLambdaUpdate TrajectoryLambdaController::observe(
    float_t lambda, float_t hpwl, float_t overflow) {
    TrajectoryLambdaUpdate update;
    update.density_step = ++density_steps_;
    update.lambda_before = lambda;
    update.lambda_after = lambda;

    if (config_.interval <= 0 || update.density_step % config_.interval != 0)
        return update;

    update.evaluated = true;
    if (!initialized_) {
        initialized_ = true;
        initial_overflow_ = std::max(overflow, config_.target_overflow);
        previous_overflow_ = overflow;
        previous_desired_ = initial_overflow_;
    }

    const float_t progress = std::min(1.0f,
        (float_t)update.density_step / std::max(1, config_.horizon_steps));
    const float_t smoothstep = progress * progress * (3.0f - 2.0f * progress);
    update.desired_overflow = config_.target_overflow +
        (initial_overflow_ - config_.target_overflow) * (1.0f - smoothstep);

    const float_t error_scale = std::max(config_.target_overflow, 1.0e-4f);
    update.controlled_overflow = overflow;
    if (config_.predictive_steps > 0.0f && progress >= 0.80f) {
        const float_t overflow_velocity = overflow - previous_overflow_;
        update.controlled_overflow = std::max(0.0f,
            overflow + config_.predictive_steps * overflow_velocity);
    }
    float_t overflow_error = update.controlled_overflow - update.desired_overflow;
    const bool in_deadband = config_.deadband > 0.0f && progress >= 1.0f &&
        std::fabs(overflow_error) <= config_.deadband;
    if (in_deadband) {
        overflow_error = 0.0f;
    } else if (config_.deadband > 0.0f && progress >= 1.0f) {
        overflow_error -= std::copysign(config_.deadband, overflow_error);
    }
    update.error = std::max(-2.0f, std::min(4.0f,
        overflow_error / error_scale));
    const float_t actual_drop = (previous_overflow_ - overflow) / error_scale;
    const float_t desired_drop =
        (previous_desired_ - update.desired_overflow) / error_scale;
    update.trend_error = std::max(-2.0f, std::min(2.0f,
        desired_drop - actual_drop));
    if (progress < 1.0f) {
        integral_error_ = std::max(-3.0f, std::min(3.0f,
            integral_error_ + update.error));
    } else {
        // Do not carry the density-ramp backlog into constrained refinement.
        // At this point the controller should regulate the final boundary,
        // not continue paying for earlier trajectory error.
        integral_error_ = 0.0f;
    }
    update.integral_error = integral_error_;

    // The trajectory supplies a moving density constraint. PI-D feedback
    // corrects position and velocity error; the small feed-forward term pays
    // for the expected increase in density force before the trajectory ends.
    const float_t kp = progress < 1.0f
        ? (update.error >= 0.0f ? config_.kp_up : config_.kp_down)
        : (update.error >= 0.0f ? config_.post_kp_up : config_.post_kp_down);
    const float_t kd = progress < 1.0f ? config_.kd : config_.post_kd;
    float_t feedforward = config_.feedforward;
    if (progress > 0.85f)
        feedforward *= std::max(0.0f, (1.0f - progress) / 0.15f);
    float_t log_delta = feedforward + kp * update.error +
        config_.ki * integral_error_ + kd * update.trend_error;
    if (in_deadband)
        log_delta = 0.0f;

    // Once density is at the final boundary, HPWL quality becomes the
    // secondary feedback signal. This is deliberately inactive while the
    // placement is still materially infeasible.
    if (progress < 1.0f && overflow <= config_.target_overflow + 0.002f) {
        const float_t hpwl_pressure = std::max(0.0f,
            hpwl / std::max(config_.hpwl_baseline, 1.0f) - 1.0f);
        log_delta -= 0.08f * hpwl_pressure;
    }

    log_delta = std::max(-config_.max_log_decrease,
                         std::min(config_.max_log_increase, log_delta));
    update.lambda_after = std::max(config_.lambda_min,
        std::min(config_.lambda_max, lambda * std::exp(log_delta)));
    if (overflow > config_.target_overflow * 1.05f)
        update.lambda_after = std::max(update.lambda_after,
                                       config_.startup_floor);

    const float_t relative_change = update.lambda_after /
        std::max(lambda, config_.lambda_min);
    if (relative_change > 1.001f)
        update.action = BaselineLambdaAction::INCREASE;
    else if (relative_change < 0.999f)
        update.action = BaselineLambdaAction::DECREASE;
    else
        update.action = BaselineLambdaAction::HOLD;

    update.hpwl_ratio = config_.hpwl_baseline /
        std::max(std::fabs(hpwl), 1.0f);
    update.overflow_ratio = config_.target_overflow /
        std::max(std::fabs(overflow), 1.0e-8f);
    previous_overflow_ = overflow;
    previous_desired_ = update.desired_overflow;
    return update;
}
