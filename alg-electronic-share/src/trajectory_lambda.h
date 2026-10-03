#pragma once

#include "baseline_lambda.h"

struct TrajectoryLambdaConfig {
    int interval = 25;
    int horizon_steps = 1150;
    float_t hpwl_baseline = 7.0e7f;
    float_t target_overflow = 0.05998f;
    float_t lambda_min = 1.0e-4f;
    float_t lambda_max = 100.0f;
    float_t startup_floor = 0.1f;
    float_t kp_up = 0.55f;
    float_t kp_down = 1.40f;
    float_t ki = 0.025f;
    float_t kd = 0.35f;
    float_t feedforward = 0.10f;
    float_t max_log_increase = 0.30f;
    float_t max_log_decrease = 0.50f;
    float_t post_kp_up = 15.0f;
    float_t post_kp_down = 8.0f;
    float_t post_kd = 0.30f;
};

struct TrajectoryLambdaUpdate {
    bool evaluated = false;
    int density_step = 0;
    float_t hpwl_ratio = 0.0f;
    float_t overflow_ratio = 0.0f;
    float_t desired_overflow = 0.0f;
    float_t error = 0.0f;
    float_t trend_error = 0.0f;
    float_t integral_error = 0.0f;
    float_t lambda_before = 0.0f;
    float_t lambda_after = 0.0f;
    BaselineLambdaAction action = BaselineLambdaAction::NONE;
};

class TrajectoryLambdaController {
public:
    explicit TrajectoryLambdaController(TrajectoryLambdaConfig config = {});

    TrajectoryLambdaUpdate observe(float_t lambda, float_t hpwl,
                                   float_t overflow);

private:
    TrajectoryLambdaConfig config_;
    int density_steps_ = 0;
    bool initialized_ = false;
    float_t initial_overflow_ = 0.0f;
    float_t previous_overflow_ = 0.0f;
    float_t previous_desired_ = 0.0f;
    float_t integral_error_ = 0.0f;
};
