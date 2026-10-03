#include "../src/baseline_lambda.h"

#include <cmath>
#include <cstdio>

namespace {

bool close(float_t lhs, float_t rhs) {
    return std::fabs(lhs - rhs) < 1.0e-7f;
}

BaselineLambdaUpdate advance(BaselineLambdaController& controller,
                             float_t lambda, float_t hpwl, float_t overflow) {
    BaselineLambdaUpdate result;
    for (int i = 0; i < 50; ++i) result = controller.observe(lambda, hpwl, overflow);
    return result;
}

}  // namespace

int main() {
    BaselineLambdaController controller;

    // HPWL is relatively better than overflow: strengthen density.
    BaselineLambdaUpdate increase = advance(controller, 0.01f, 7.0e7f, 0.12f);
    if (!increase.evaluated || increase.action != BaselineLambdaAction::INCREASE ||
        !close(increase.lambda_after, 0.02f)) {
        std::fprintf(stderr, "expected lambda increase\n");
        return 1;
    }

    // Overflow is relatively better than HPWL: recover wirelength.
    BaselineLambdaUpdate decrease = advance(controller, 0.02f, 1.4e8f, 0.06f);
    if (!decrease.evaluated || decrease.action != BaselineLambdaAction::DECREASE ||
        !close(decrease.lambda_after, 0.01f)) {
        std::fprintf(stderr, "expected lambda decrease\n");
        return 2;
    }

    // Similar normalized quality lies inside the deadband.
    BaselineLambdaUpdate hold = advance(controller, 0.01f, 7.0e7f, 0.06f);
    if (!hold.evaluated || hold.action != BaselineLambdaAction::HOLD ||
        !close(hold.lambda_after, 0.01f)) {
        std::fprintf(stderr, "expected lambda hold\n");
        return 3;
    }

    BaselineLambdaConfig guard_cfg;
    guard_cfg.interval = 1;
    guard_cfg.guarded = true;
    guard_cfg.recovery_overflow = 0.057f;
    BaselineLambdaController guard(guard_cfg);

    BaselineLambdaUpdate spread = guard.observe(0.5f, 7.0e7f, 0.20f);
    if (spread.action != BaselineLambdaAction::INCREASE ||
        !close(spread.lambda_after, 2.0f)) return 4;
    BaselineLambdaUpdate recover = guard.observe(2.0f, 1.2e8f, 0.04f);
    if (recover.action != BaselineLambdaAction::DECREASE ||
        !(recover.lambda_after < 2.0f && recover.lambda_after >= 1.0f)) return 5;
    BaselineLambdaController balanced_guard(guard_cfg);
    BaselineLambdaUpdate guard_hold = balanced_guard.observe(
        recover.lambda_after, 1.0e8f, 0.057f);
    if (guard_hold.action == BaselineLambdaAction::INCREASE) return 6;
    BaselineLambdaUpdate respread = guard.observe(1.0e-4f, 9.0e7f, 0.061f);
    if (respread.action != BaselineLambdaAction::INCREASE ||
        !(respread.lambda_after > 1.0e-4f)) return 7;

    std::printf("baseline lambda controller test passed\n");
    return 0;
}
