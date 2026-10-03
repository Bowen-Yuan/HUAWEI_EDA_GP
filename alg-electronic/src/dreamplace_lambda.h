#pragma once

#include "baseline_lambda.h"

struct DreamplaceLambdaConfig {
    int interval = 1;
    float_t stop_overflow = 0.10f;
    float_t ref_hpwl = 350000.0f;
    float_t lower_pcof = 0.95f;
    float_t upper_pcof = 1.05f;
    float_t lambda_min = 1.0e-4f;
    float_t lambda_max = 100.0f;
};

struct DreamplaceLambdaUpdate {
    bool evaluated = false;
    bool frozen = false;
    int density_step = 0;
    float_t delta_hpwl = 0.0f;
    float_t normalized_delta = 0.0f;
    float_t multiplier = 1.0f;
    float_t lambda_before = 0.0f;
    float_t lambda_after = 0.0f;
    BaselineLambdaAction action = BaselineLambdaAction::NONE;
};

class DreamplaceLambdaController {
public:
    explicit DreamplaceLambdaController(DreamplaceLambdaConfig config = {});

    DreamplaceLambdaUpdate observe(float_t lambda, float_t hpwl,
                                   float_t overflow);

private:
    DreamplaceLambdaConfig config_;
    int density_steps_ = 0;
    bool initialized_ = false;
    bool frozen_ = false;
    float_t previous_hpwl_ = 0.0f;
};
