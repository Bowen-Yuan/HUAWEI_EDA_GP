#pragma once

#include "epsilon_active/batch_acceptance.hpp"
#include "epsilon_active/bisection.hpp"
#include "epsilon_active/coarse_flow.hpp"
#include "epsilon_active/compact_recovery.hpp"
#include "epsilon_active/density_coordinate.hpp"
#include "epsilon_active/lambda_controller.hpp"
#include "epsilon_active/optimizer.hpp"
#include "epsilon_active/recovery.hpp"
#include "epsilon_active/swap_recovery.hpp"
#include "epsilon_active/transport.hpp"
#include "epsilon_active/types.hpp"

#include <cstdint>
#include <filesystem>
#include <vector>

namespace ea {

struct NetBatchConfig {
    // Zero preserves the existing node-wise direction exactly. A positive
    // value adds a net-coherent search direction; it is never an objective
    // term and every trial is still audited with exact costs.
    Real weight = 0.0;
    int degree_limit = 64;
};

// Exact overlap is still evaluated node-by-node.  This optional search
// heuristic shares the measured overlap direction across the movable pins of
// each incident net, which can activate nodes whose individual exact
// subgradient is zero.  It changes only the trial direction, never the
// objective or acceptance oracle.
struct DensityNetShareConfig {
    Real weight = 0.0;
    int degree_limit = 64;
    Real zero_gradient_threshold = 1.0e-12;
    // Remove the component of the shared density direction parallel to the
    // current exact HPWL subgradient.  This is a direction preconditioner,
    // not a penalty or a change to either objective.
    Real hpwl_orthogonal_projection = 0.0;
    int until_iteration = -1;
};

struct RegionalPriceConfig {
    bool enabled = false;
    int update_interval = 10;
    Real rho = 0.25;
    Real target_excess = 0.0;
    Real min_price = 0.25;
    Real max_price = 4.0;
    std::vector<Real> initial_prices;
};

// Optional in-process continuation switch.  The objective and all exact
// oracles remain unchanged; this only changes the first-order step scale and
// density multiplier after the selected iteration.  A negative iteration
// keeps the historical single-stage behavior exactly.
struct LateStageConfig {
    int switch_iteration = -1;
    Real step_multiplier = 1.0;
    Real lambda_multiplier = 1.0;
    bool reset_optimizer = false;
};

struct PlaceConfig {
    int bins_x = 512;
    int bins_y = 512;
    int iterations = 1000;
    int threads = 0;
    int degree_limit = 100;
    int log_every = 10;
    Real target_density = 1.0;
    Real stop_overflow = 0.07;
    Real hpwl_epsilon = 125.0;
    Real density_epsilon = 0.0;
    Real active_power = 4.0;
    Real density_active_power = 1.0;
    Real step_fraction = 0.003;
    Real max_step_multiplier = 4.0;
    BatchAcceptanceConfig batch_acceptance;
    OptimizerKind optimizer = OptimizerKind::Adam;
    Real beta1 = 0.90;
    Real beta2 = 0.99;
    Real momentum = 0.70;
    Real numerical_epsilon = 1.0e-8;
    LambdaConfig lambda;
    bool adaptive_epsilon = true;
    int adaptive_interval = 25;
    int adaptive_window = 50;
    Real adaptive_gain = 0.5;
    Real adaptive_deadband = 0.0025;
    Real adaptive_min_scale = 0.70;
    Real adaptive_max_scale = 1.10;
    Real adaptive_max_log_step = 0.03;
    Real lower_overflow = 0.065;
    NetBatchConfig net_batch;
    DensityNetShareConfig density_net_share;
    RegionalPriceConfig regional_price;
    LateStageConfig late_stage;
    BisectionConfig bisection;
    CoarseFlowConfig coarse_flow;
    TransportConfig transport;
    RecoveryConfig recovery;
    DensityCoordinateConfig density_coordinate;
    CompactRecoveryConfig compact_recovery;
    SwapRecoveryConfig swap_recovery;
    RecoveryConfig post_recovery;
    SwapRecoveryConfig post_swap_recovery;
    SwapRecoveryConfig assignment_recovery;
    int snapshot_every = 0;
    std::filesystem::path output_dir;
};

struct PlaceResult {
    IterationMetrics final_metrics;
    bool feasible = false;
    int selected_iteration = -1;
    int objective_evaluations = 0;
    int accepted_batches = 0;
    int rejected_batches = 0;
    int batch_trials = 0;
    BisectionStats bisection;
    CoarseFlowStats coarse_flow;
    TransportStats transport;
    RecoveryStats recovery;
    DensityCoordinateStats density_coordinate;
    CompactRecoveryStats compact_recovery;
    SwapRecoveryStats swap_recovery;
    RecoveryStats post_recovery;
    SwapRecoveryStats post_swap_recovery;
    SwapRecoveryStats assignment_recovery;
    double wall_seconds = 0.0;
    std::vector<Real> final_regional_prices;
};

PlaceResult global_place(Database& db, const PlaceConfig& config);

}  // namespace ea
