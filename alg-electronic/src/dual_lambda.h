#pragma once

#include "baseline_lambda.h"

struct DualLambdaConfig {
    int interval = 25;
    float_t target_overflow = 0.0995f;
    float_t step_size = 1.0f;
    float_t lambda_min = 1.0e-4f;
    float_t lambda_max = 100.0f;
};

struct DualLambdaUpdate {
    bool evaluated = false;
    int density_step = 0;
    float_t constraint_violation = 0.0f;
    float_t lambda_before = 0.0f;
    float_t lambda_after = 0.0f;
    BaselineLambdaAction action = BaselineLambdaAction::NONE;
};

class DualLambdaController {
public:
    explicit DualLambdaController(DualLambdaConfig config = {});

    DualLambdaUpdate observe(float_t lambda, float_t overflow);

private:
    DualLambdaConfig config_;
    int density_steps_ = 0;
};

