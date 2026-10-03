#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <numeric>

// ============================================================================
// Bayesian Optimization for Lambda Selection
//
// Uses a Gaussian Process (GP) with Matern 5/2 kernel to model the objective
// landscape f(λ) = HPWL_norm + penalty * [max(0, overflow - target)]².
//
// At each decision point, Expected Improvement (EI) selects the next λ.
// This naturally balances exploration and exploitation.
// ============================================================================

struct LambdaBOConfig {
    double length_scale   = 1.5;
    double signal_std     = 0.5;
    double noise_std      = 0.05;
    double penalty_weight = 10.0;
    int    window_size    = 40;
    int    refit_interval = 5;
    double lambda_min     = 1e-4;
    double lambda_max     = 50.0;
    int    grid_points    = 200;
    double ei_xi          = 0.01;
};

class LambdaBayesianOpt {
public:
    LambdaBayesianOpt();
    explicit LambdaBayesianOpt(const LambdaBOConfig& cfg);

    // Core BO interface: call once per solver iteration
    double step(double current_lambda, double hpwl, double overflow,
                double target_overflow);

    // Force-select a new λ (called at phase transitions)
    double force_select(double current_lambda, double hpwl, double overflow,
                        double target_overflow);

    void reset();

    bool   is_active()       const { return n_obs_ >= 5; }
    double best_lambda()     const { return best_lambda_; }
    double best_objective()  const { return best_obj_; }

private:
    LambdaBOConfig cfg_;

    std::vector<double> X_;       // log(λ)
    std::vector<double> y_;       // objective values
    int n_obs_ = 0;
    int iter_since_refit_ = 0;

    std::vector<double> alpha_;   // K⁻¹y
    std::vector<double> L_;       // Cholesky factor (flattened)
    bool gp_valid_ = false;

    double best_lambda_ = 0.1;
    double best_obj_ = 1e30;

    double kernel(double x1, double x2) const;
    void fit_gp();
    void predict(double x_star, double& mean, double& std) const;
    double expected_improvement(double x, double y_best) const;
    double optimize_ei(double y_best) const;
    double compute_objective(double hpwl, double overflow, double target) const;
};
