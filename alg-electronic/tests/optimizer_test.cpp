#include "../src/optimizer.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "optimizer_test failed: %s\n", message);
        std::exit(1);
    }
}

float_t step(OptimizerKind kind, float_t gradient, int age,
             const OptimizerHyperparameters& hp,
             float_t& m, float_t& v, float_t& vmax) {
    return optimizer_coordinate_step(kind, gradient, 2.0f, age,
                                     hp, m, v, vmax);
}

}  // namespace

int main() {
    OptimizerHyperparameters hp;
    hp.momentum = 0.70f;
    hp.beta1 = 0.90f;
    hp.beta2 = 0.99f;
    hp.step_scale = 1.0f;
    hp.max_step_multiplier = 10.0f;

    float_t m = 0.0f, v = 0.0f, vmax = 0.0f;
    const float_t hb = step(OptimizerKind::HeavyBall, 1.0f, 1, hp, m, v, vmax);
    require(std::fabs(hb - 0.6f) < 1.0e-5f,
            "heavy-ball must preserve the legacy EMA update");

    m = v = vmax = 0.0f;
    const float_t nag = step(OptimizerKind::Nesterov, 1.0f, 1, hp, m, v, vmax);
    require(nag > hb && nag < 2.0f,
            "Nesterov must react more strongly than heavy-ball on a fresh gradient");

    m = v = vmax = 0.0f;
    const float_t adam = step(OptimizerKind::Adam, 2.0f, 1, hp, m, v, vmax);
    require(std::fabs(adam - 2.0f) < 1.0e-4f,
            "bias-corrected Adam first step must have unit normalized scale");

    m = v = vmax = 0.0f;
    (void)step(OptimizerKind::AMSGrad, 4.0f, 1, hp, m, v, vmax);
    const float_t old_vmax = vmax;
    (void)step(OptimizerKind::AMSGrad, 0.5f, 2, hp, m, v, vmax);
    require(vmax >= old_vmax, "AMSGrad denominator must be monotone");

    m = v = vmax = 0.0f;
    const float_t rms = step(OptimizerKind::RMSProp, -3.0f, 1, hp, m, v, vmax);
    require(rms < 0.0f, "RMSProp must preserve descent direction sign");

    m = v = vmax = 0.0f;
    const float_t ada1 = std::fabs(step(OptimizerKind::AdaGrad, 1.0f, 1, hp, m, v, vmax));
    const float_t ada2 = std::fabs(step(OptimizerKind::AdaGrad, 1.0f, 2, hp, m, v, vmax));
    require(ada2 < ada1, "AdaGrad coordinate step must decay after accumulation");

    std::puts("optimizer_test passed");
    return 0;
}
