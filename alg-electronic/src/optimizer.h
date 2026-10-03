#pragma once

#include "types.h"

struct OptimizerHyperparameters {
    float_t momentum = 0.70f;
    float_t beta1 = 0.90f;
    float_t beta2 = 0.99f;
    float_t epsilon = 1.0e-6f;
    float_t step_scale = 1.0f;
    float_t max_step_multiplier = 4.0f;
};

float_t optimizer_coordinate_step(
    OptimizerKind kind,
    float_t gradient,
    float_t learning_rate,
    int age,
    const OptimizerHyperparameters& hp,
    float_t& first_moment,
    float_t& second_moment,
    float_t& max_second_moment);

bool optimizer_uses_adaptive_state(OptimizerKind kind);
const char* optimizer_name(OptimizerKind kind);
