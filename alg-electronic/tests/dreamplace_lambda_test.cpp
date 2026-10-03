#include "../src/dreamplace_lambda.h"

#include <cstdio>
#include <cstdlib>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "dreamplace_lambda_test failed: %s\n", message);
        std::exit(1);
    }
}

}  // namespace

int main() {
    DreamplaceLambdaConfig cfg;
    DreamplaceLambdaController controller(cfg);

    auto first = controller.observe(0.02f, 90.0e6f, 0.30f);
    require(first.evaluated && first.action == BaselineLambdaAction::HOLD,
            "first observation must establish the HPWL reference");

    auto improving = controller.observe(0.02f, 89.8e6f, 0.25f);
    require(improving.multiplier > 1.0f && improving.lambda_after > 0.02f,
            "decreasing HPWL above stop overflow must increase density weight");

    auto worsening = controller.observe(improving.lambda_after, 91.0e6f, 0.20f);
    require(worsening.multiplier >= cfg.lower_pcof &&
            worsening.multiplier <= cfg.upper_pcof,
            "paper multiplier must stay inside the RePlAce bounds");

    auto stop = controller.observe(worsening.lambda_after, 90.9e6f, 0.099f);
    require(stop.frozen && stop.lambda_after == worsening.lambda_after,
            "controller must freeze at the official overflow threshold");

    auto frozen = controller.observe(stop.lambda_after, 80.0e6f, 0.20f);
    require(frozen.frozen && !frozen.evaluated &&
            frozen.lambda_after == stop.lambda_after,
            "DREAMPlace update mask must not reactivate after freezing");

    std::puts("dreamplace_lambda_test passed");
    return 0;
}
