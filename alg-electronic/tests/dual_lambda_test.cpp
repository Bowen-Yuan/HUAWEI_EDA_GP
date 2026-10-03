#include "../src/dual_lambda.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>

static void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "dual_lambda_test failed: %s\n", message);
        std::exit(1);
    }
}

int main() {
    DualLambdaConfig cfg;
    cfg.interval = 2;
    cfg.target_overflow = 0.10f;
    cfg.step_size = 2.0f;
    cfg.lambda_min = 0.01f;
    cfg.lambda_max = 1.0f;
    DualLambdaController controller(cfg);

    auto skipped = controller.observe(0.20f, 0.30f);
    require(!skipped.evaluated, "controller evaluated before its interval");

    auto raised = controller.observe(0.20f, 0.30f);
    require(raised.evaluated, "controller skipped its update interval");
    require(raised.action == BaselineLambdaAction::INCREASE,
            "positive violation did not increase lambda");
    require(std::fabs(raised.lambda_after - 0.60f) < 1.0e-6f,
            "dual increase has the wrong magnitude");

    controller.observe(raised.lambda_after, 0.05f);
    auto lowered = controller.observe(raised.lambda_after, 0.05f);
    require(lowered.action == BaselineLambdaAction::DECREASE,
            "negative violation did not decrease lambda");
    require(std::fabs(lowered.lambda_after - 0.50f) < 1.0e-6f,
            "dual decrease has the wrong magnitude");

    DualLambdaConfig clipped_cfg = cfg;
    clipped_cfg.interval = 1;
    clipped_cfg.step_size = 100.0f;
    DualLambdaController clipped(clipped_cfg);
    require(std::fabs(clipped.observe(0.20f, 1.0f).lambda_after - 1.0f) < 1.0e-6f,
            "upper projection failed");
    require(std::fabs(clipped.observe(0.20f, 0.0f).lambda_after - 0.01f) < 1.0e-6f,
            "lower projection failed");

    std::puts("dual_lambda_test passed");
    return 0;
}
