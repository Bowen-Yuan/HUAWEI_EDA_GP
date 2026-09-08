#pragma once

#include "epsilon_active/density.hpp"
#include "epsilon_active/optimizer.hpp"

#include <functional>
#include <string>
#include <vector>

namespace ea {

// Shared mechanics only: individual modules own how the auxiliary direction
// is constructed.  No candidate is accepted, rejected, or restored here.
struct NonlocalDescentConfig {
    int iterations = 100;
    Real learning_rate = 0.002;
    Real beta1 = 0.90;
    Real beta2 = 0.999;
    Real numerical_epsilon = 1.0e-8;
    Real maximum_delta = 4.0;
    OptimizerKind optimizer = OptimizerKind::AMSGrad;
    Real hpwl_epsilon = 0.5;
    Real exact_density_weight = 1.0;
    Real auxiliary_density_weight = 1.0;
    Real lambda_initial = 0.25;
    Real lambda_minimum = 0.01;
    Real lambda_maximum = 256.0;
    Real start_overflow_percent = 15.0;
    Real stop_overflow_percent = 7.0;
    int lambda_update_interval = 5;
    Real lambda_kp = 0.20;
    Real lambda_ki = 0.02;
    Real lambda_kd = 0.10;
};

struct NonlocalDescentStats {
    int objective_evaluations = 0;
    int nonzero_auxiliary_steps = 0;
    Real maximum_overflow = 0.0;
};

using AuxiliaryDirection = std::function<void(
    const Database&, ExactOverlapDensity&, std::vector<Real>&,
    std::vector<Real>&)>;
using NonlocalTelemetry = std::function<void(int, Real, const DensityMetrics&, Real,
                                             Real, Real, Real)>;

NonlocalDescentStats run_nonlocal_descent(Database& db, int bins_x, int bins_y,
    Real target_density, const NonlocalDescentConfig& config,
    const AuxiliaryDirection& auxiliary, const NonlocalTelemetry& telemetry = {});

}  // namespace ea
