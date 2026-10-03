#include "optimizer.h"

#include <algorithm>
#include <cmath>

namespace {

float_t clamp_step(float_t step, float_t learning_rate,
                   const OptimizerHyperparameters& hp) {
    const float_t cap = hp.max_step_multiplier * learning_rate;
    return std::max(-cap, std::min(cap, step));
}

}  // namespace

bool optimizer_uses_adaptive_state(OptimizerKind kind) {
    return kind != OptimizerKind::HeavyBall && kind != OptimizerKind::Nesterov;
}

const char* optimizer_name(OptimizerKind kind) {
    switch (kind) {
    case OptimizerKind::HeavyBall: return "heavy-ball";
    case OptimizerKind::Adam: return "adam";
    case OptimizerKind::Nesterov: return "nesterov-ema";
    case OptimizerKind::RMSProp: return "rmsprop";
    case OptimizerKind::AMSGrad: return "amsgrad";
    case OptimizerKind::AdaGrad: return "adagrad";
    }
    return "unknown";
}

float_t optimizer_coordinate_step(
    OptimizerKind kind,
    float_t gradient,
    float_t learning_rate,
    int age,
    const OptimizerHyperparameters& hp,
    float_t& first_moment,
    float_t& second_moment,
    float_t& max_second_moment) {
    const int safe_age = std::max(1, age);

    if (kind == OptimizerKind::HeavyBall) {
        first_moment = hp.momentum * first_moment +
                       (1.0f - hp.momentum) * gradient;
        return learning_rate * first_moment;
    }

    if (kind == OptimizerKind::Nesterov) {
        first_moment = hp.momentum * first_moment +
                       (1.0f - hp.momentum) * gradient;
        const float_t direction = hp.momentum * first_moment +
                                  (1.0f - hp.momentum) * gradient;
        return learning_rate * direction;
    }

    if (kind == OptimizerKind::AdaGrad) {
        second_moment += gradient * gradient;
        const float_t step = hp.step_scale * learning_rate * gradient /
            (std::sqrt(second_moment) + hp.epsilon);
        return clamp_step(step, learning_rate, hp);
    }

    second_moment = hp.beta2 * second_moment +
                    (1.0f - hp.beta2) * gradient * gradient;
    const float_t b2_correction = std::max(
        1.0e-8f, 1.0f - std::pow(hp.beta2, (float_t)safe_age));
    const float_t corrected_second = second_moment / b2_correction;

    if (kind == OptimizerKind::RMSProp) {
        const float_t step = hp.step_scale * learning_rate * gradient /
            (std::sqrt(corrected_second) + hp.epsilon);
        return clamp_step(step, learning_rate, hp);
    }

    first_moment = hp.beta1 * first_moment +
                   (1.0f - hp.beta1) * gradient;
    const float_t b1_correction = std::max(
        1.0e-8f, 1.0f - std::pow(hp.beta1, (float_t)safe_age));
    const float_t corrected_first = first_moment / b1_correction;

    float_t denominator_second = corrected_second;
    if (kind == OptimizerKind::AMSGrad) {
        max_second_moment = std::max(max_second_moment, corrected_second);
        denominator_second = max_second_moment;
    }
    const float_t step = hp.step_scale * learning_rate * corrected_first /
        (std::sqrt(denominator_second) + hp.epsilon);
    return clamp_step(step, learning_rate, hp);
}
