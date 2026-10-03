#pragma once

#include "types.h"

enum class BaselineLambdaAction {
    NONE,
    INCREASE,
    DECREASE,
    HOLD,
};

struct BaselineLambdaConfig {
    int interval = 50;
    float_t hpwl_baseline = 7.0e7f;
    float_t overflow_baseline = 0.06f;
    float_t increase_factor = 2.0f;
    float_t decrease_factor = 0.5f;
    float_t hysteresis = 0.05f;
    float_t lambda_min = 1.0e-4f;
    float_t lambda_max = 20.0f;
    bool guarded = false;
    float_t recovery_overflow = 0.0595f;
    float_t spread_lambda_floor = 2.0f;
};

struct BaselineLambdaUpdate {
    bool evaluated = false;
    int density_step = 0;
    float_t hpwl_ratio = 0.0f;
    float_t overflow_ratio = 0.0f;
    float_t lambda_before = 0.0f;
    float_t lambda_after = 0.0f;
    BaselineLambdaAction action = BaselineLambdaAction::NONE;
};

class BaselineLambdaController {
public:
    explicit BaselineLambdaController(BaselineLambdaConfig config = {});

    BaselineLambdaUpdate observe(float_t lambda, float_t hpwl, float_t overflow);

private:
    BaselineLambdaConfig config_;
    int density_steps_ = 0;
    float_t overflow_error_ema_ = 0.0f;
};

const char* baseline_lambda_action_name(BaselineLambdaAction action);
