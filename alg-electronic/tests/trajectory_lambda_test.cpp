#include "../src/trajectory_lambda.h"

#include <cstdio>

int main() {
    TrajectoryLambdaConfig cfg;
    cfg.interval = 1;
    cfg.horizon_steps = 10;
    TrajectoryLambdaController controller(cfg);

    TrajectoryLambdaUpdate first = controller.observe(0.01f, 7.0e7f, 0.50f);
    if (!first.evaluated || first.lambda_after < cfg.startup_floor ||
        first.desired_overflow >= 0.50f) return 1;

    TrajectoryLambdaUpdate behind = controller.observe(
        first.lambda_after, 7.5e7f, 0.50f);
    if (behind.action != BaselineLambdaAction::INCREASE ||
        behind.lambda_after <= first.lambda_after) return 2;

    TrajectoryLambdaController fast_controller(cfg);
    TrajectoryLambdaUpdate seed = fast_controller.observe(1.0f, 8.0e7f, 0.50f);
    TrajectoryLambdaUpdate ahead = fast_controller.observe(
        seed.lambda_after, 8.0e7f, 0.20f);
    if (ahead.action != BaselineLambdaAction::DECREASE ||
        ahead.lambda_after >= seed.lambda_after) return 3;

    TrajectoryLambdaConfig boundary_cfg;
    boundary_cfg.interval = 1;
    boundary_cfg.horizon_steps = 1;
    TrajectoryLambdaController boundary(boundary_cfg);
    TrajectoryLambdaUpdate feasible = boundary.observe(
        10.0f, 8.0e7f, boundary_cfg.target_overflow - 0.001f);
    if (feasible.action != BaselineLambdaAction::DECREASE) return 4;
    TrajectoryLambdaUpdate infeasible = boundary.observe(
        feasible.lambda_after, 8.0e7f, boundary_cfg.target_overflow + 0.001f);
    if (infeasible.action != BaselineLambdaAction::INCREASE) return 5;

    TrajectoryLambdaConfig predictive_cfg;
    predictive_cfg.interval = 1;
    predictive_cfg.horizon_steps = 2;
    predictive_cfg.predictive_steps = 1.0f;
    predictive_cfg.deadband = 0.001f;
    TrajectoryLambdaController predictive(predictive_cfg);
    TrajectoryLambdaUpdate predictive_seed = predictive.observe(
        1.0f, 8.0e7f, predictive_cfg.target_overflow + 0.020f);
    TrajectoryLambdaUpdate predictive_brake = predictive.observe(
        predictive_seed.lambda_after, 8.0e7f,
        predictive_cfg.target_overflow + 0.010f);
    if (predictive_brake.controlled_overflow >
            predictive_cfg.target_overflow + 1.0e-6f ||
        predictive_brake.action != BaselineLambdaAction::HOLD) return 6;

    TrajectoryLambdaConfig deadband_cfg;
    deadband_cfg.interval = 1;
    deadband_cfg.horizon_steps = 1;
    deadband_cfg.deadband = 0.001f;
    TrajectoryLambdaController deadband(deadband_cfg);
    TrajectoryLambdaUpdate held = deadband.observe(
        2.0f, 8.0e7f, deadband_cfg.target_overflow + 0.0005f);
    if (held.action != BaselineLambdaAction::HOLD ||
        std::fabs(held.lambda_after - 2.0f) > 1.0e-6f) return 7;

    std::printf("trajectory lambda controller test passed\n");
    return 0;
}
