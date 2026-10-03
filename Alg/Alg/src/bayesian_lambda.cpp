#include "bayesian_lambda.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <cassert>

// ============================================================================
// Bayesian Optimization for Lambda Selection — Implementation
// ============================================================================

LambdaBayesianOpt::LambdaBayesianOpt() : cfg_(LambdaBOConfig{}) {
    reset();
}

LambdaBayesianOpt::LambdaBayesianOpt(const LambdaBOConfig& cfg) : cfg_(cfg) {
    reset();
}

// Matern 5/2 kernel: smooth enough to model HPWL/overflow landscape,
// but not infinitely differentiable (more realistic for optimization landscapes)
// k(r) = σ_f² · (1 + √5·r/ρ + 5r²/3ρ²) · exp(-√5·r/ρ)
double LambdaBayesianOpt::kernel(double x1, double x2) const {
    double r = std::abs(x1 - x2) / cfg_.length_scale;
    double sqrt5 = 2.23606797749979;
    double tr = sqrt5 * r;
    double k = cfg_.signal_std * cfg_.signal_std
             * (1.0 + tr + tr * tr / 3.0) * std::exp(-tr);
    return k;
}

// Compute objective: we want to MINIMIZE this
// f(λ) = HPWL/1e7 + w · [max(0, overflow - target)]²
// Normalization by 1e7 keeps values in a reasonable range for GP
double LambdaBayesianOpt::compute_objective(double hpwl, double overflow,
                                             double target) const {
    double over = std::max(0.0, overflow - target);
    return hpwl / 1.0e7 + cfg_.penalty_weight * over * over;
}

// ---------------------------------------------------------------------------
// GP inference via Cholesky decomposition
// K = L·Lᵀ, α = K⁻¹y = L⁻ᵀ(L⁻¹y)
// ---------------------------------------------------------------------------
void LambdaBayesianOpt::fit_gp() {
    int n = n_obs_;
    if (n < 2) { gp_valid_ = false; return; }

    // Build kernel matrix K + σ_n²·I
    std::vector<double> K(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double kval = kernel(X_[i], X_[j]);
            if (i == j) kval += cfg_.noise_std * cfg_.noise_std;
            K[i * n + j] = kval;
            K[j * n + i] = kval;  // symmetric
        }
    }

    // Cholesky decomposition K = L·Lᵀ
    L_.resize(n * n);
    std::fill(L_.begin(), L_.end(), 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = K[i * n + j];
            for (int k = 0; k < j; ++k)
                sum -= L_[i * n + k] * L_[j * n + k];

            if (i == j) {
                if (sum <= 1e-12) {
                    // Numerical issue — add jitter
                    sum = 1e-10;
                }
                L_[i * n + i] = std::sqrt(sum);
            } else {
                L_[i * n + j] = sum / L_[j * n + j];
            }
        }
    }

    // Solve L·α_tmp = y  (forward substitution)
    alpha_.resize(n);
    for (int i = 0; i < n; ++i) {
        double sum = y_[i];
        for (int j = 0; j < i; ++j)
            sum -= L_[i * n + j] * alpha_[j];
        alpha_[i] = sum / L_[i * n + i];
    }

    // Solve Lᵀ·α = α_tmp  (backward substitution)
    for (int i = n - 1; i >= 0; --i) {
        double sum = alpha_[i];
        for (int j = i + 1; j < n; ++j)
            sum -= L_[j * n + i] * alpha_[j];
        alpha_[i] = sum / L_[i * n + i];
    }

    gp_valid_ = true;
}

// Predict mean and std at x_star
void LambdaBayesianOpt::predict(double x_star, double& mean, double& std_dev) const {
    if (!gp_valid_ || n_obs_ < 2) {
        mean = 0.0;
        std_dev = cfg_.signal_std;
        return;
    }

    int n = n_obs_;

    // k_* = K(x_star, X)
    std::vector<double> k_star(n);
    for (int i = 0; i < n; ++i)
        k_star[i] = kernel(x_star, X_[i]);

    // k_** = K(x_star, x_star)
    double k_ss = kernel(x_star, x_star) + cfg_.noise_std * cfg_.noise_std;

    // mean = k_* · α
    mean = 0.0;
    for (int i = 0; i < n; ++i)
        mean += k_star[i] * alpha_[i];

    // v = L⁻¹ · k_*
    std::vector<double> v(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = k_star[i];
        for (int j = 0; j < i; ++j)
            sum -= L_[i * n + j] * v[j];
        v[i] = sum / L_[i * n + i];
    }

    // σ² = k_** - vᵀv
    double var = k_ss;
    for (int i = 0; i < n; ++i)
        var -= v[i] * v[i];

    if (var < 1e-12) var = 1e-12;
    std_dev = std::sqrt(var);
}

// ---------------------------------------------------------------------------
// Expected Improvement (for minimization)
// EI(x) = (y_best - μ(x)) · Φ(z) + σ(x) · φ(z)
// where z = (y_best - μ(x) - ξ) / σ(x)
//
// ξ (cfg_.ei_xi) controls exploration: ξ>0 encourages more exploration
// ---------------------------------------------------------------------------
double LambdaBayesianOpt::expected_improvement(double x, double y_best) const {
    double mu, sigma;
    predict(x, mu, sigma);

    double diff = y_best - mu - cfg_.ei_xi;
    if (sigma < 1e-12) return (diff > 0) ? diff : 0.0;

    double z = diff / sigma;

    // Standard normal CDF and PDF (using error function)
    double phi = 0.5 * (1.0 + std::erf(z / 1.4142135623730951));  // Φ(z)
    double pdf = std::exp(-0.5 * z * z) / 2.5066282746310002;      // φ(z)

    return diff * phi + sigma * pdf;
}

// Grid-search for max EI in log-λ space
double LambdaBayesianOpt::optimize_ei(double y_best) const {
    double best_x = best_lambda_;
    double best_ei = -1e30;

    double lo = std::log(cfg_.lambda_min);
    double hi = std::log(cfg_.lambda_max);

    for (int i = 0; i < cfg_.grid_points; ++i) {
        double frac = (double)i / (cfg_.grid_points - 1);
        double x = lo + frac * (hi - lo);
        double ei = expected_improvement(x, y_best);

        if (ei > best_ei) {
            best_ei = ei;
            best_x = x;
        }
    }

    return std::exp(best_x);
}

// ---------------------------------------------------------------------------
// Main BO step
// ---------------------------------------------------------------------------
double LambdaBayesianOpt::step(double current_lambda, double hpwl,
                                double overflow, double target_overflow) {
    double obj = compute_objective(hpwl, overflow, target_overflow);

    // Add observation to sliding window
    double lx = std::log(std::max(cfg_.lambda_min, current_lambda));
    if (n_obs_ >= cfg_.window_size) {
        // Remove oldest
        X_.erase(X_.begin());
        y_.erase(y_.begin());
        n_obs_--;
    }
    X_.push_back(lx);
    y_.push_back(obj);
    n_obs_++;

    // Track best
    if (obj < best_obj_) {
        best_obj_ = obj;
        best_lambda_ = current_lambda;
    }

    iter_since_refit_++;

    // Only refit and select new λ periodically
    if (iter_since_refit_ < cfg_.refit_interval || n_obs_ < 5) {
        return current_lambda;  // keep current λ
    }

    iter_since_refit_ = 0;

    // Fit GP
    fit_gp();
    if (!gp_valid_) return current_lambda;

    // Optimize EI to get next λ
    double new_lambda = optimize_ei(best_obj_);

    // Trust region: don't change λ too drastically in one step
    double ratio = new_lambda / current_lambda;
    if (ratio > 3.0)  new_lambda = current_lambda * 3.0;
    if (ratio < 0.33) new_lambda = current_lambda * 0.33;

    // Clamp
    new_lambda = std::max(cfg_.lambda_min, std::min(cfg_.lambda_max, new_lambda));

    return new_lambda;
}

// Force a new λ selection (for phase transitions)
double LambdaBayesianOpt::force_select(double current_lambda, double hpwl,
                                        double overflow, double target_overflow) {
    // Record the current state first
    double obj = compute_objective(hpwl, overflow, target_overflow);
    double lx = std::log(std::max(cfg_.lambda_min, current_lambda));

    if (n_obs_ >= cfg_.window_size) {
        X_.erase(X_.begin());
        y_.erase(y_.begin());
        n_obs_--;
    }
    X_.push_back(lx);
    y_.push_back(obj);
    n_obs_++;

    if (obj < best_obj_) {
        best_obj_ = obj;
        best_lambda_ = current_lambda;
    }

    iter_since_refit_ = 0;

    // Force GP refit
    fit_gp();
    if (!gp_valid_ || n_obs_ < 5) {
        // Not enough data — perturb current λ slightly
        return std::max(cfg_.lambda_min,
                        std::min(cfg_.lambda_max, current_lambda * 1.5));
    }

    double new_lambda = optimize_ei(best_obj_);
    double ratio = new_lambda / current_lambda;
    if (ratio > 3.0)  new_lambda = current_lambda * 3.0;
    if (ratio < 0.33) new_lambda = current_lambda * 0.33;
    new_lambda = std::max(cfg_.lambda_min, std::min(cfg_.lambda_max, new_lambda));

    return new_lambda;
}

void LambdaBayesianOpt::reset() {
    X_.clear();
    y_.clear();
    alpha_.clear();
    L_.clear();
    n_obs_ = 0;
    iter_since_refit_ = 0;
    gp_valid_ = false;
    best_lambda_ = 0.1;
    best_obj_ = 1e30;
}
