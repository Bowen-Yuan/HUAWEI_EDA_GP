#include "dreamplace_lambda.h"

#include <algorithm>
#include <cmath>

DreamplaceLambdaController::DreamplaceLambdaController(
    DreamplaceLambdaConfig config) : config_(config) {}

DreamplaceLambdaUpdate DreamplaceLambdaController::observe(
    float_t lambda, float_t hpwl, float_t overflow) {
    DreamplaceLambdaUpdate update;
    update.density_step = ++density_steps_;
    update.lambda_before = lambda;
    update.lambda_after = lambda;

    if (frozen_) {
        update.frozen = true;
        return update;
    }
    if (config_.interval <= 0 || update.density_step % config_.interval != 0)
        return update;

    update.evaluated = true;
    if (overflow <= config_.stop_overflow) {
        frozen_ = true;
        update.frozen = true;
        update.action = BaselineLambdaAction::HOLD;
        return update;
    }

    if (!initialized_) {
        initialized_ = true;
        previous_hpwl_ = hpwl;
        update.action = BaselineLambdaAction::HOLD;
        return update;
    }

    update.delta_hpwl = hpwl - previous_hpwl_;
    update.normalized_delta = update.delta_hpwl /
        std::max(config_.ref_hpwl, 1.0f);
    if (update.normalized_delta < 0.0f) {
        update.multiplier = config_.upper_pcof * std::max(
            std::pow(0.9999f, (float_t)update.density_step), 0.98f);
    } else {
        update.multiplier = std::pow(
            config_.upper_pcof, 1.0f - update.normalized_delta);
        update.multiplier = std::max(config_.lower_pcof,
            std::min(config_.upper_pcof, update.multiplier));
    }

    update.lambda_after = std::max(config_.lambda_min,
        std::min(config_.lambda_max, lambda * update.multiplier));
    if (update.lambda_after > lambda * 1.001f)
        update.action = BaselineLambdaAction::INCREASE;
    else if (update.lambda_after < lambda * 0.999f)
        update.action = BaselineLambdaAction::DECREASE;
    else
        update.action = BaselineLambdaAction::HOLD;
    previous_hpwl_ = hpwl;
    return update;
}
