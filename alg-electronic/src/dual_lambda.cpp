#include "dual_lambda.h"

#include <algorithm>

DualLambdaController::DualLambdaController(DualLambdaConfig config)
    : config_(config) {}

DualLambdaUpdate DualLambdaController::observe(
    float_t lambda, float_t overflow) {
    DualLambdaUpdate update;
    update.density_step = ++density_steps_;
    update.lambda_before = lambda;
    update.lambda_after = lambda;

    if (config_.interval <= 0 || update.density_step % config_.interval != 0)
        return update;

    update.evaluated = true;
    update.constraint_violation = overflow - config_.target_overflow;
    update.lambda_after = std::max(config_.lambda_min,
        std::min(config_.lambda_max,
                 lambda + config_.step_size * update.constraint_violation));

    const float_t delta = update.lambda_after - lambda;
    if (delta > 1.0e-6f)
        update.action = BaselineLambdaAction::INCREASE;
    else if (delta < -1.0e-6f)
        update.action = BaselineLambdaAction::DECREASE;
    else
        update.action = BaselineLambdaAction::HOLD;
    return update;
}

