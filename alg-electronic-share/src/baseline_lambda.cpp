#include "baseline_lambda.h"

#include <algorithm>
#include <cmath>

BaselineLambdaController::BaselineLambdaController(BaselineLambdaConfig config)
    : config_(config) {}

BaselineLambdaUpdate BaselineLambdaController::observe(
    float_t lambda, float_t hpwl, float_t overflow) {
    BaselineLambdaUpdate update;
    update.density_step = ++density_steps_;
    update.lambda_before = lambda;
    update.lambda_after = lambda;

    if (config_.interval <= 0 || update.density_step % config_.interval != 0)
        return update;

    update.evaluated = true;
    const float_t safe_hpwl = std::max(std::fabs(hpwl), 1.0f);
    const float_t safe_overflow = std::max(std::fabs(overflow), 1.0e-8f);
    update.hpwl_ratio = config_.hpwl_baseline / safe_hpwl;
    update.overflow_ratio = config_.overflow_baseline / safe_overflow;

    if (config_.guarded) {
        // Smooth primal-dual feedback around a slightly conservative overflow
        // target.  The previous bang-bang guard repeatedly dropped lambda to
        // its minimum and then jumped to the spread floor; those large pulses
        // destroyed HPWL progress.  This controller keeps the long-term
        // overflow error in an EMA and changes lambda continuously.
        const float_t target = std::max(1.0e-6f, config_.recovery_overflow);
        const float_t overflow_error = (overflow - target) /
                                       std::max(config_.overflow_baseline, 1.0e-6f);
        overflow_error_ema_ = 0.75f * overflow_error_ema_ +
                              0.25f * overflow_error;
        const float_t hpwl_pressure = std::max(0.0f,
            hpwl / std::max(config_.hpwl_baseline, 1.0f) - 1.0f);

        // Density error raises lambda; excessive HPWL gently lowers it when
        // density is already feasible. Limit each update to [0.5x, 2x].
        float_t log_factor = 1.35f * overflow_error +
                             0.35f * overflow_error_ema_;
        if (overflow <= config_.overflow_baseline)
            log_factor -= 0.12f * hpwl_pressure;
        log_factor = std::max(std::log(std::max(config_.decrease_factor, 1.0e-3f)),
                     std::min(std::log(std::max(config_.increase_factor, 1.001f)),
                              log_factor));
        update.lambda_after = std::max(config_.lambda_min,
            std::min(config_.lambda_max, lambda * std::exp(log_factor)));

        // If a recovery pulse has reduced lambda nearly to zero and overflow
        // is materially infeasible, restart from a modest active value.
        if (overflow > config_.overflow_baseline * 1.05f)
            update.lambda_after = std::max(update.lambda_after,
                                           config_.spread_lambda_floor);

        const float_t relative_change = update.lambda_after /
                                        std::max(lambda, config_.lambda_min);
        if (relative_change > 1.001f)
            update.action = BaselineLambdaAction::INCREASE;
        else if (relative_change < 0.999f)
            update.action = BaselineLambdaAction::DECREASE;
        else
            update.action = BaselineLambdaAction::HOLD;
        return update;
    }

    // A larger baseline/current ratio means that metric is relatively closer
    // to (or better than) its baseline. Spend weight on the other objective.
    if (update.hpwl_ratio > update.overflow_ratio * (1.0f + config_.hysteresis)) {
        update.lambda_after = std::min(config_.lambda_max,
            lambda * config_.increase_factor);
        update.action = BaselineLambdaAction::INCREASE;
    } else if (update.overflow_ratio >
               update.hpwl_ratio * (1.0f + config_.hysteresis)) {
        update.lambda_after = std::max(config_.lambda_min,
            lambda * config_.decrease_factor);
        update.action = BaselineLambdaAction::DECREASE;
    } else {
        update.action = BaselineLambdaAction::HOLD;
    }
    return update;
}

const char* baseline_lambda_action_name(BaselineLambdaAction action) {
    switch (action) {
    case BaselineLambdaAction::INCREASE: return "increase";
    case BaselineLambdaAction::DECREASE: return "decrease";
    case BaselineLambdaAction::HOLD: return "hold";
    default: return "none";
    }
}
