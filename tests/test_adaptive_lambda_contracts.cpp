// Synthetic contract tests for the adaptive_lambda_gp search policy.
// No benchmark data is required; every scenario is a tiny in-memory layout.

#include "adaptive_lambda_gp.hpp"

#include "epsilon_active/density.hpp"
#include "epsilon_active/hpwl.hpp"
#include "epsilon_active/optimizer.hpp"
#include "epsilon_active/step_policy.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nsgp::adaptive;
using ea::Real;

namespace {

int g_failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++g_failures;
        std::cerr << "FAILED: " << message << '\n';
    }
}

void check_close(Real a, Real b, Real tolerance, const std::string& message) {
    check(std::abs(a - b) <= tolerance, message);
}

ea::Database three_cell_row() {
    // Three 10x10 cells on a 30x10 die with three 10x10 bins.  Density is
    // exactly 1.0 everywhere, so any HPWL-converging move must create excess
    // in the middle bin while lowering HPWL.
    ea::Database db;
    db.xl = 0.0;
    db.yl = 0.0;
    db.xh = 30.0;
    db.yh = 10.0;
    db.movable_area = 300.0;
    db.nodes = {{0, "a", 5.0, 5.0, 10.0, 10.0, false, false, "N"},
                {1, "b", 15.0, 5.0, 10.0, 10.0, false, false, "N"},
                {2, "c", 25.0, 5.0, 10.0, 10.0, false, false, "N"}};
    db.movable_ids = {0, 1, 2};
    db.node_pin_count = {1, 1, 1};
    db.node_pin_offsets = {0, 1, 2, 3};
    db.node_pin_indices = {0, 1, 2};
    db.pins = {{0, 0.0, 0.0}, {1, 0.0, 0.0}, {2, 0.0, 0.0}};
    db.nets = {{0, "n", 0, 3, 1.0}};
    return db;
}

AdaptiveLambdaOptions restore_options() {
    AdaptiveLambdaOptions options;
    // Five iterations keep the first proposal inside the phase-A corridor;
    // later funnel phases reject everything and force the fallback.
    options.iterations = 5;
    options.active_ensemble = false;
    options.optimizer = "sgd";
    options.step_policy = "constant";
    options.learning_rate = 1.0;
    options.maximum_delta_bins = 10000.0;
    options.lambda.fixed = true;
    options.lambda.initial = 0.25;
    return options;
}

void test_step_policy_contracts() {
    ea::StepControllerConfig config;

    // Constant: identical decisions every iteration, bins converted to units.
    {
        auto controller = ea::make_step_controller(ea::parse_step_policy("constant"), config);
        controller->reset(10, 2.0);
        for (int iteration = 0; iteration <= 10; ++iteration) {
            const auto decision = controller->propose(iteration);
            check_close(decision.learning_rate, 0.02, 0.0, "constant lr invariant");
            check_close(decision.maximum_delta, 1.0, 1.0e-12, "constant bins conversion");
        }
    }

    // Cosine: lr(0) == lr0, monotone non-increasing, lr(end) == lr_min.
    {
        auto controller = ea::make_step_controller(ea::parse_step_policy("cosine"), config);
        controller->reset(10, 1.0);
        Real previous = 0.0;
        for (int iteration = 0; iteration <= 10; ++iteration) {
            const auto decision = controller->propose(iteration);
            if (iteration == 0) {
                check_close(decision.learning_rate, 0.02, 1.0e-12, "cosine lr(0) == lr0");
            }
            if (iteration > 0) {
                check(decision.learning_rate <= previous + 1.0e-15, "cosine monotone");
            }
            previous = decision.learning_rate;
        }
        check_close(previous, 0.002, 1.0e-9, "cosine lr(end) == lr_min");
    }

    // Trust: growth after consecutive accepts, shrink on reject, clamped.
    {
        auto controller = ea::make_step_controller(ea::parse_step_policy("trust"), config);
        controller->reset(100, 1.0);
        check_close(controller->propose(0).maximum_delta, 1.0, 0.0, "trust initial radius");
        for (int i = 0; i < 3; ++i) controller->observe({true, 0});
        check_close(controller->propose(1).maximum_delta, 1.25, 1.0e-12, "trust grows after 3 accepts");
        controller->observe({false, 9});
        check_close(controller->propose(2).maximum_delta, 0.625, 1.0e-12, "trust shrinks on reject");
        for (int i = 0; i < 40; ++i) controller->observe({false, 9});
        check_close(controller->propose(3).maximum_delta,
                    config.radius_min_bins, 0.0, "trust radius clamped at minimum");
        controller->reset(100, 1.0);
        for (int i = 0; i < 40; ++i) controller->observe({true, 0});
        check_close(controller->propose(4).maximum_delta,
                    config.radius_max_bins, 0.0, "trust radius clamped at maximum");
    }

    // reset_streaks clears the growth streak without changing the radius.
    {
        auto controller = ea::make_step_controller(ea::parse_step_policy("trust"), config);
        controller->reset(100, 1.0);
        controller->observe({true, 0});
        controller->observe({true, 0});
        controller->reset_streaks();
        controller->observe({true, 0});
        controller->observe({true, 0});
        check_close(controller->propose(0).maximum_delta, 1.0, 1.0e-12,
                    "reset_streaks defers growth");
        controller->observe({true, 0});
        check_close(controller->propose(0).maximum_delta, 1.25, 1.0e-12,
                    "growth after restored streak");
    }
}

void test_lambda_controller_contracts() {
    FunnelLambdaConfig defaults;  // 15% explore, 7% final, hold 0.20, end 0.75

    // Case A: far below the phase-A corridor must not increase lambda.
    {
        FunnelLambdaController controller(defaults);
        const Real ratio = controller.maybe_update(5, 50, 0.072);
        check(ratio < 1.0, "case A: below corridor lowers lambda");
    }

    // Case B: above a 7% corridor must increase lambda.
    {
        FunnelLambdaConfig config = defaults;
        config.explore_overflow_percent = 7.0;  // constant 7% corridor
        FunnelLambdaController controller(config);
        const Real ratio = controller.maybe_update(5, 50, 0.085);
        check(ratio > 1.0, "case B: above corridor raises lambda");
    }

    // Case C: with gains isolated to the trend term, a worsening trend must
    // raise lambda more than an improving trend at the same error.
    {
        FunnelLambdaConfig config = defaults;
        config.explore_overflow_percent = 7.0;
        config.kp = 0.0;
        config.ki = 0.0;
        config.kd = 0.10;
        FunnelLambdaController worsening(config);
        worsening.maybe_update(5, 50, 0.080);
        const Real lambda_worse = [&] {
            worsening.maybe_update(10, 50, 0.085);
            return worsening.lambda();
        }();
        FunnelLambdaController improving(config);
        improving.maybe_update(5, 50, 0.090);
        const Real lambda_better = [&] {
            improving.maybe_update(10, 50, 0.085);
            return improving.lambda();
        }();
        check(lambda_worse > lambda_better,
              "case C: worsening trend outruns improving trend");
    }

    // Case D: final-lock boost applies only at p >= contract_end.
    {
        FunnelLambdaConfig config = defaults;
        config.explore_overflow_percent = 7.0;  // constant 7% corridor
        FunnelLambdaController mid(config);
        mid.maybe_update(50, 100, 0.08);
        FunnelLambdaController locked(config);
        locked.maybe_update(80, 100, 0.08);
        check(locked.lambda() > mid.lambda(), "case D: final-lock boost");
        check_close(locked.lambda() / mid.lambda(), std::exp(0.10), 1.0e-9,
                    "case D: boost magnitude is exp(extra)");
    }

    // Case E: extreme sustained errors saturate inside [min, max].
    {
        FunnelLambdaConfig config = defaults;
        config.explore_overflow_percent = 7.0;
        FunnelLambdaController controller(config);
        for (int iteration = config.update_interval;
             iteration <= 40 * config.update_interval;
             iteration += config.update_interval) {
            controller.maybe_update(iteration, 2000, 0.50);
            check(controller.lambda() >= config.minimum - 1.0e-15 &&
                      controller.lambda() <= config.maximum + 1.0e-15,
                  "case E: lambda within bounds");
        }
        check_close(controller.lambda(), config.maximum, 1.0e-12,
                    "case E: saturates at maximum");
        for (int iteration = 41 * config.update_interval;
             iteration <= 120 * config.update_interval;
             iteration += config.update_interval) {
            controller.maybe_update(iteration, 2000, 0.0001);
        }
        check_close(controller.lambda(), config.minimum, 1.0e-12,
                    "case E: recovers to minimum");
    }
}

void test_acceptance_contracts() {
    FunnelAcceptanceConfig config;

    check(funnel_accept(config, 0.09, 100.0, {100.0, 0.07}, {99.0, 0.085}) ==
              FunnelDecision::ExploreAccept,
          "funnel: explore accepts hpwl gain inside corridor");
    check(funnel_accept(config, 0.08, 100.0, {100.0, 0.085}, {100.03, 0.082}) ==
              FunnelDecision::RecoveryAccept,
          "funnel: recovery accepts bounded overflow repair");
    check(funnel_accept(config, 0.15, 100.0, {100.0, 0.07}, {99.0, 0.17}) ==
              FunnelDecision::Reject,
          "funnel: hard ceiling rejects");
    check(funnel_accept(config, 0.08, 100.0, {100.0, 0.085}, {100.0, 0.09}) ==
              FunnelDecision::Reject,
          "funnel: worsening overflow outside corridor rejects");
    check(funnel_accept(config, 0.09, 100.0, {100.0, 0.07},
                        {std::numeric_limits<Real>::quiet_NaN(), 0.08}) ==
              FunnelDecision::Reject,
          "funnel: non-finite candidate rejects");
    check(funnel_accept(config, 0.09, 100.0, {100.0, 0.07}, {100.0, 0.08}) ==
              FunnelDecision::Reject,
          "funnel: hpwl must strictly decrease in explore mode");
}

void test_best_feasible_restore_contract() {
    // One iteration accepts an exploring candidate at ~11.5% overflow; the
    // stage must still return the iteration-0 feasible layout.
    {
        ea::Database db = three_cell_row();
        const auto options = restore_options();
        const auto stats = run_adaptive_lambda_search(db, 3, 1, 1.0, options, false);
        check(stats.accepted >= 1, "restore: exploring candidate accepted");
        check(stats.last_before_restore.overflow > 0.07,
              "restore: search ended above the final cap");
        check(stats.last_before_restore.overflow <= 0.16,
              "restore: excursion stayed under the hard ceiling");
        check_close(stats.best_feasible.hpwl, 20.0, 1.0e-9,
                    "restore: best feasible stays at the input");
        check_close(stats.final_selected.hpwl, 20.0, 1.0e-9,
                    "restore: stage output restores input hpwl");
        check_close(stats.final_selected.overflow, 0.0, 1.0e-15,
                    "restore: stage output restores input overflow");
        check_close(db.nodes[0].x, 5.0, 1.0e-12, "restore: node a position");
        check_close(db.nodes[2].x, 25.0, 1.0e-12, "restore: node c position");
    }

    // If the input itself is infeasible under the final target the stage
    // must fail loudly instead of emitting a >target layout.
    {
        ea::Database db = three_cell_row();
        auto options = restore_options();
        options.lambda.final_overflow_percent = 0.1;
        bool threw = false;
        try {
            run_adaptive_lambda_search(db, 3, 1, 0.995, options, false);
        } catch (const std::runtime_error&) {
            threw = true;
        }
        check(threw, "restore: infeasible input under final target throws");
    }
}

void test_optimizer_step_smoke() {
    const std::vector<std::string> optimizers = {
        "adam", "amsgrad", "adagrad", "heavy-ball", "sgd", "normalized-sgd",
        "dual-averaging"};
    const std::vector<std::string> steps = {"constant", "cosine", "trust"};
    const std::vector<Real> gradient = {1.0, -1.0, 0.5, 0.25};

    for (const std::string& optimizer_name : optimizers) {
        for (const std::string& step_name : steps) {
            auto make = [&](std::vector<std::vector<Real>>& history) {
                auto optimizer = ea::make_optimizer(
                    ea::parse_optimizer(optimizer_name), 0.9, 0.999, 0.9, 1.0e-8);
                optimizer->reset(gradient.size());
                ea::StepControllerConfig config;
                config.learning_rate = 0.02;
                auto step = ea::make_step_controller(
                    ea::parse_step_policy(step_name), config);
                step->reset(10, 1.0);
                std::vector<Real> delta;
                for (int iteration = 1; iteration <= 10; ++iteration) {
                    const auto decision = step->propose(iteration);
                    optimizer->compute_delta(gradient, decision.learning_rate,
                                             decision.maximum_delta, delta);
                    for (const Real value : delta) {
                        check(std::isfinite(value),
                              optimizer_name + "+" + step_name + " finite delta");
                        check(std::abs(value) <= decision.maximum_delta * (1.0 + 1.0e-12),
                              optimizer_name + "+" + step_name + " bounded delta");
                    }
                    history.push_back(delta);
                    step->observe({iteration % 3 != 0, 0});
                }
            };
            std::vector<std::vector<Real>> first_pass;
            make(first_pass);
            std::vector<std::vector<Real>> second_pass;
            make(second_pass);
            check(first_pass == second_pass,
                  optimizer_name + "+" + step_name + " reset reproducibility");
        }
    }
}

void test_validation_rejects_invalid_options() {
    AdaptiveLambdaOptions options;
    auto expect_throw = [&](const char* label) {
        bool threw = false;
        try {
            options.validate();
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, std::string("validation rejects ") + label);
    };
    options.iterations = 0;
    expect_throw("iterations < 1");
    options = restore_options();
    options.lambda.explore_overflow_percent = 6.0;
    expect_throw("explore below final");
    options = restore_options();
    options.lambda.hard_overflow_percent = 6.9;
    expect_throw("hard below explore");
    options = restore_options();
    options.lambda.hold_fraction = 0.8;
    expect_throw("hold beyond contract end");
    options = restore_options();
    options.optimizer = "bogus";
    expect_throw("unknown optimizer");
    options = restore_options();
    options.step_policy = "bogus";
    expect_throw("unknown step policy");
    bool threw = false;
    try {
        (void)ea::parse_step_policy("bogus");
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "parse_step_policy rejects unknown names");
}

}  // namespace

int main() {
    test_step_policy_contracts();
    test_lambda_controller_contracts();
    test_acceptance_contracts();
    test_best_feasible_restore_contract();
    test_optimizer_step_smoke();
    test_validation_rejects_invalid_options();
    if (g_failures) {
        std::cerr << g_failures << " adaptive lambda contract checks failed\n";
        return 1;
    }
    std::cout << "adaptive lambda contracts passed\n";
    return 0;
}
