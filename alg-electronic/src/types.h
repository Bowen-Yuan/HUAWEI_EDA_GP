#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <cmath>
#include <cfloat>
#include <algorithm>

using float_t = float;

struct Cell {
    int id;
    float_t x, y;         // current position
    float_t width, height;
    bool is_terminal;     // true = fixed IO pad
    float_t area;
};

struct Net {
    int id;
    std::vector<int> cell_ids;  // indices into cells array
    float_t hpwl_x;
    float_t hpwl_y;
};

struct BinGrid {
    int nx, ny;
    float_t bin_w, bin_h;
    float_t chip_xl, chip_yl, chip_xh, chip_yh;
    float_t target_density;
    std::vector<float_t> density;     // nx * ny
    std::vector<float_t> overflow;    // nx * ny
    float_t total_overflow;
    float_t available_area = 0.0f;
};

struct PlacementResult {
    float_t hpwl;
    float_t overflow;
    int iterations;
    double runtime_sec;
};

struct VisualizationConfig {
    bool enabled = false;
    int every = 0;
    std::string output_dir;
};

enum class OptimizerKind {
    HeavyBall,
    Adam,
    Nesterov,
    RMSProp,
    AMSGrad,
    AdaGrad
};

struct SolverConfig {
    int max_iters;           // Max outer iterations
    float_t step_scale;      // Base step size scaling
    float_t step_min;        // Minimum step size
    float_t step_max;        // Maximum step size
    float_t momentum;        // Heavy-ball momentum coefficient
    float_t lambda_init;     // Initial density penalty weight
    float_t lambda_growth;   // λ growth factor when overflow is high
    float_t lambda_decay;    // λ decay factor when overflow is low
    float_t target_overflow; // Target overflow ratio
    float_t grad_clip;       // Gradient clipping threshold
    int p0_iters;            // Phase 0: HPWL-only iterations (default 150)
    int coarse_iters;        // Phase 1: iterations with coarse bins (default 100)
    int medium_iters;        // Phase 2: iterations with medium bins (default 150)
    float_t p1_lambda_init;  // Phase 1: initial λ for density penalty (default 0.005)
    bool use_bayesian_opt;   // Use Bayesian optimization for λ selection
    int lambda_update_interval = 50;
    float_t hpwl_baseline = 7.0e7f;
    float_t overflow_baseline = 0.06f;
    float_t overflow_limit = 0.06f;
    float_t lambda_increase_factor = 2.0f;
    float_t lambda_decrease_factor = 0.5f;
    float_t lambda_ratio_hysteresis = 0.05f;
    float_t lambda_control_min = 1.0e-4f;
    float_t lambda_control_max = 20.0f;
    float_t target_density_override = -1.0f;
    bool use_adam = false;
    OptimizerKind optimizer_kind = OptimizerKind::HeavyBall;
    float_t adam_beta1 = 0.80f;
    float_t adam_beta2 = 0.99f;
    float_t adam_epsilon = 1.0e-6f;
    float_t adam_step_scale = 0.70f;
    float_t adam_max_step_multiplier = 4.0f;
    int adam_late_step = -1;
    float_t adam_late_step_scale = -1.0f;
    int amsgrad_switch_to_adam_step = -1;
    bool use_lambda_guard = true;
    bool use_trajectory_lambda = false;
    bool use_dual_lambda = false;
    bool use_dreamplace_lambda = false;
    float_t dual_lambda_step_size = 1.0f;
    float_t dreamplace_lambda_ref_hpwl = 350000.0f;
    int dreamplace_lambda_interval = 1;
    float_t lambda_guard_low = 0.0598f;
    float_t lambda_guard_floor = 0.1f;
    int lambda_guard_interval = 25;
    int lambda_trajectory_horizon = 1150;
    float_t lambda_trajectory_target = 0.05998f;
    float_t lambda_trajectory_predictive_steps = 0.0f;
    float_t lambda_trajectory_deadband = 0.0f;
    float_t fine_momentum = 0.70f;
    float_t fine_cooldown_floor = 0.15f;
    bool use_projected_hpwl = false;
    bool use_adaptive_restart = false;
    float_t projection_overflow_margin = 0.002f;
    bool use_hpwl_bundle = false;
    int hpwl_bundle_size = 4;
    int hpwl_bundle_interval = 5;
    bool hpwl_bundle_adaptive_cuts = false;
    int hpwl_bundle_min_interval = 2;
    float_t hpwl_bundle_cosine_trigger = 0.90f;
    float_t hpwl_bundle_prox_scale = 1.0f;
    float_t hpwl_bundle_current_mix = 0.35f;
    bool hpwl_bundle_extrema_mix = false;
    float_t hpwl_bundle_mix_gain = 0.30f;
    int hpwl_bundle_mix_decay_steps = 0;
    float_t hpwl_bundle_start_overflow = 1.0f;
    bool lambda_scale_current_hpwl = false;
    bool gradient_diagnostics = false;
    bool use_filter_trust = false;
    int filter_interval = 5;
    float_t filter_hpwl_tolerance = 2.0e-4f;
    float_t filter_overflow_tolerance = 2.0e-4f;
    float_t filter_reject_scale = 0.50f;
    float_t filter_growth = 1.05f;
    float_t filter_min_step_scale = 0.125f;
    bool use_feasible_trust = false;
    float_t feasible_trust_slack = 1.0e-4f;
    float_t feasible_trust_hpwl_tolerance = 0.0f;
    float_t feasible_trust_reject_scale = 0.5f;
    float_t feasible_trust_growth = 1.02f;
    float_t feasible_trust_min_step_scale = 0.01f;
    float_t merit_trust_hpwl_target = -1.0f;
    float_t merit_trust_start_overflow = 1.0f;
    float_t merit_trust_tolerance = 0.0f;
    float_t merit_trust_reject_scale = 0.5f;
    float_t merit_trust_growth = 1.02f;
    float_t merit_trust_min_step_scale = 0.001f;
    VisualizationConfig visualization;
};

struct CWTRConfig {
    int max_sweeps = 15;
    int candidates_per_dim = 5;
    float_t trust_radius_initial = 2.0f;
    float_t overflow_relax = 1.05f;
    float_t hpwl_improve_min = 0.9999f;
    int swap_attempts = 500;
    float_t target_density = 0.80f;
    float_t overflow_limit = 0.06f;
};
