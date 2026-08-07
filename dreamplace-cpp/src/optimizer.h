#pragma once

#include "electric.h"
#include "types.h"

#include <cstdint>
#include <string>
#include <vector>

namespace dpcpp {

enum class WirelengthModel { WeightedAverage, ExactHpwl };
enum class GlobalOptimizer { DreamplaceNesterov, HeavyBall, Adam, AMSGrad, AdaGrad };
enum class LambdaPolicy { Dreamplace, Trajectory, Ratio };

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
    int bundle_interval = 5;
    int bundle_min_interval = 2;
    Real bundle_cosine_trigger = 0.90;
    Real bundle_start_overflow = 0.12;
    Real bundle_prox_scale = 1.0;
    Real bundle_current_mix = 0.50;
    bool feasible_refinement = false;
    Real refinement_lower_overflow = 0.065;
    Real refinement_learning_rate_scale = 0.05;
    Real refinement_lambda_gain = 0.20;
    GlobalOptimizer refinement_optimizer = GlobalOptimizer::AMSGrad;
    int snapshot_every = 0;
    std::string snapshot_dir;
    std::string metrics_path;
};

struct GlobalPlaceResult {
    Metrics final_metrics;
    Metrics best_feasible_metrics;
    bool have_feasible = false;
    int objective_evaluations = 0;
    double wall_time_seconds = 0.0;
};

GlobalPlaceResult global_place(Database& db, std::vector<Filler>& fillers,
                               const GlobalPlaceConfig& config);

}  // namespace dpcpp
