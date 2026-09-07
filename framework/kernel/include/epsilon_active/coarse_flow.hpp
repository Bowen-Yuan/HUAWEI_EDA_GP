#pragma once

#include "epsilon_active/density.hpp"
#include "epsilon_active/types.hpp"

namespace ea {

struct CoarseFlowConfig {
    bool enabled = false;
    int bins_x = 64;
    int bins_y = 64;
    int passes = 1;
    bool exact_anchors = false;
    // Zero keeps the exhaustive legacy behavior. Positive values bound the
    // global candidate generator for runtime-safe screening.
    int max_source_bins = 0;
    int max_nodes_per_source = 0;
    int anchor_grid_stride = 1;
    bool connected_commodities = false;
    // Try several fractions of the coarse rigid displacement for each
    // connected commodity. Disabled by default to preserve H348 behavior.
    bool commodity_multiscale = false;
    // Also screen x-only and y-only rigid directions for each flow edge.
    bool commodity_axis_split = false;
    bool commodity_hpwl_guard = false;
    // Jointly select disjoint commodity proposals and audit the whole batch.
    bool joint_assignment = false;
    bool joint_allow_nonmonotone = false;
    int joint_proposals_per_edge = 8;
    bool joint_ignore_edge_area = false;
    Real joint_hpwl_budget_ratio = 0.0;
    int commodity_degree_limit = 16;
    int commodity_max_nodes = 8;
};

struct CoarseFlowStats {
    int passes = 0;
    int flow_edges = 0;
    int moves = 0;
    int exact_anchor_trials = 0;
    int exact_anchor_rejections = 0;
    Real planned_area = 0.0;
    Real moved_area = 0.0;
    Real initial_coarse_overflow = 0.0;
    Real final_coarse_overflow = 0.0;
    Real initial_hpwl = 0.0;
    Real initial_overflow = 0.0;
    Real final_hpwl = 0.0;
    Real final_overflow = 0.0;
    double wall_seconds = 0.0;
    int commodity_groups = 0;
    int commodity_moves = 0;
    int commodity_scale_trials = 0;
    int joint_batches = 0;
    int joint_batch_trials = 0;
};

CoarseFlowStats coarse_capacity_flow(
    Database& db, ExactOverlapDensity& density, const CoarseFlowConfig& config);

}  // namespace ea
