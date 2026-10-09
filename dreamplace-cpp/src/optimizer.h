#pragma once

#include "electric.h"
#include "types.h"

#include <cstdint>
#include <string>
#include <vector>

namespace dpcpp {

enum class WirelengthModel { WeightedAverage, ExactHpwl };
enum class GlobalOptimizer { DreamplaceNesterov, HeavyBall, Adam, AMSGrad, AdaGrad };
enum class LambdaPolicy { Dreamplace, Trajectory, Ratio, BandDual };

const char* wirelength_model_name(WirelengthModel model);
const char* global_optimizer_name(GlobalOptimizer optimizer);
const char* lambda_policy_name(LambdaPolicy policy);

struct GlobalPlaceConfig {
    int bins_x = 512;
    int bins_y = 512;
    int iterations = 1000;
    int degree_limit = 100;
    int log_every = 10;
    Real target_density = 1.0;
    Real stop_overflow = 0.07;
    Real density_weight_scale = 8.0e-5;
    Real gamma_scale = 4.0;
    Real initial_learning_rate = 0.01;
    Real gp_noise_ratio = 0.025;
    Real reference_hpwl_delta = 3.5e5;
    std::uint64_t seed = 1000;
    bool enable_fillers = true;
    WirelengthModel wirelength_model = WirelengthModel::WeightedAverage;
    GlobalOptimizer optimizer = GlobalOptimizer::DreamplaceNesterov;
    LambdaPolicy lambda_policy = LambdaPolicy::Dreamplace;
    int hpwl_only_iterations = 0;
    Real nonsmooth_step_fraction = 0.002;
    Real nonsmooth_hpwl_step_fraction = 0.005;
    Real momentum = 0.70;
    Real beta1 = 0.90;
    Real beta2 = 0.99;
    Real epsilon = 1.0e-8;
    Real max_step_multiplier = 4.0;
    int lambda_update_interval = 25;
    int lambda_trajectory_horizon = 550;
    Real hpwl_baseline = 0.0;
    Real overflow_baseline = 0.07;
    Real lambda_control_min = 1.0e-6;
    Real lambda_control_max = 1.0e12;
    bool enable_bundle = false;
    int bundle_size = 4;
    int bundle_groups = 1;
    int bundle_group_size = 3;
    int bundle_interval = 5;
    int bundle_min_interval = 2;
    Real bundle_cosine_trigger = 0.90;
    Real bundle_start_overflow = 0.12;
    Real bundle_prox_scale = 1.0;
    Real bundle_current_mix = 0.50;
    bool feasible_refinement = false;
    Real refinement_start_overflow = -1.0;
    Real refinement_lower_overflow = 0.065;
    Real refinement_learning_rate_scale = 0.05;
    Real refinement_lambda_gain = 0.20;
    GlobalOptimizer refinement_optimizer = GlobalOptimizer::AMSGrad;
    Real refinement_active_set_radius = -1.0;
    int refinement_active_set_decay_iterations = 0;
    bool tangent_refinement = false;
    int refinement_filter_backtracks = 6;
    bool legal_checkpoint_selection = false;
    int legal_checkpoint_interval = 20;
    int legal_checkpoint_detailed_passes = 1;
    bool late_legal_projection = false;
    int legal_projection_interval = 40;
    Real legal_projection_mix = 0.25;
    Real legal_projection_max_hpwl_ratio = 1.08;
    Real legal_row_force = 0.0;
    bool mixed_spectral_field = false;
    bool multilevel_density = false;
    int multilevel_min_bins = 128;
    Real multilevel_middle_overflow = 0.30;
    Real multilevel_fine_overflow = 0.15;
    int gradient_sampling_samples = 0;
    int gradient_sampling_interval = 5;
    Real gradient_sampling_radius = 0.25;
    int refinement_gradient_sampling_samples = 0;
    Real refinement_gradient_sampling_radius = 0.25;
    Real active_set_radius = 0.0;
    Real active_set_power = 2.0;
    // Experimental, opt-in radius controller.  The legacy fixed-radius path
    // remains unchanged when this flag is false.
    bool adaptive_active_set = false;
    bool adaptive_active_smart = false;
    // Experimental second-generation controller.  It keeps the exact HPWL
    // objective and only changes the epsilon-active trial direction.
    bool adaptive_active_predictive = false;
    Real adaptive_active_set_min_scale = 0.40;
    Real adaptive_active_set_max_scale = 1.20;
    int adaptive_active_set_interval = 25;
    int adaptive_active_set_window = 50;
    Real adaptive_active_set_gain = 1.0;
    Real adaptive_active_set_deadband = 0.0015;
    Real adaptive_active_set_max_log_step = 0.06;
    bool adaptive_active_set_refinement = false;
    Real adaptive_active_set_refinement_max_log_step = 0.01;
    Real adaptive_active_set_span_cap = 0.0;
    Real primal_dual_step = 0.0;
    bool serious_bundle = false;
    Real serious_step_ratio = 0.10;
    Real bundle_overflow_tolerance = 0.0;
    bool progressive_legalization = false;
    int progressive_density_iterations = 400;
    Real progressive_overflow_lower = 0.07;
    Real progressive_overflow_upper = 0.08;
    Real progressive_hpwl_target = 7.5e7;
    Real progressive_dual_kp = 0.055;
    Real progressive_dual_ki = 0.002;
    Real progressive_hpwl_gain = 0.035;
    Real progressive_dual_switch_overflow = 0.11;
    Real progressive_obstacle_scale = 0.20;
    Real progressive_row_force = 0.03;
    Real progressive_segment_force = 0.05;
    Real progressive_congestion_gain = 2.0;
    int progressive_obstacle_iterations = 200;
    int progressive_filter_backtracks = 6;
    bool progressive_filter = true;
    int snapshot_every = 0;
    std::string snapshot_dir;
    std::string metrics_path;
    bool profile = false;
};

struct GlobalPlaceResult {
    Metrics final_metrics;
    Metrics best_feasible_metrics;
    bool have_feasible = false;
    int objective_evaluations = 0;
    double wall_time_seconds = 0.0;
    Real selected_legal_hpwl = 0.0;
    int selected_legal_iteration = -1;
};

GlobalPlaceResult global_place(Database& db, std::vector<Filler>& fillers,
                               const GlobalPlaceConfig& config);

}  // namespace dpcpp
