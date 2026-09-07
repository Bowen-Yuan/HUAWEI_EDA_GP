#include "epsilon_active/bookshelf.hpp"
#include "epsilon_active/placer.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

namespace {

struct Options {
    std::filesystem::path benchmark;
    std::filesystem::path output;
    ea::PlaceConfig place;
    std::filesystem::path initial_placement;
    std::uint64_t seed = 1000;
    ea::Real sigma_ratio = 0.001;
    bool multilevel = false;
    int level_iterations = 200;
};

std::string value_after(int& index, int argc, char** argv) {
    if (++index >= argc) throw std::invalid_argument("missing command-line value");
    return argv[index];
}

void print_help() {
    std::cout
        << "epsilon_active --benchmark BASE --output DIR [options]\n"
        << "  --initial-placement FILE.pl       (optional GP checkpoint seed)\n"
        << "  --bins N | --bins-x N --bins-y N\n"
        << "  --multilevel --level-iterations N  (64/128/256/512 exact stages)\n"
        << "  --iterations N --threads N          (0=runtime default, maximum 40)\n"
        << "  --optimizer adam|amsgrad|adagrad|heavy-ball|sgd\n"
        << "  --lambda-policy dreamplace|trajectory|ratio\n"
        << "  --hpwl-epsilon X --density-epsilon X --active-power X\n"
        << "  --density-active-power X\n"
        << "  --step-fraction X --density-weight-scale X\n"
        << "  --late-switch-iteration N --late-step-multiplier X\n"
        << "  --late-lambda-multiplier X --late-reset-optimizer\n"
        << "  --net-batch-weight X --net-batch-degree-limit N\n"
        << "  --density-net-share X --density-net-degree-limit N\n"
        << "  --density-net-zero-threshold X\n"
        << "  --density-hpwl-orthogonal-projection X\n"
        << "  --density-net-share-until N\n"
        << "  --regional-price --regional-price-interval N --regional-price-rho X\n"
        << "  --regional-price-min X --regional-price-max X\n"
        << "  --batch-backtracking --batch-infeasible --batch-hpwl-budget X\n"
        << "  --batch-trials N --batch-shrink X\n"
        << "  --target-density X --stop-overflow X --degree-limit N\n"
        << "  --adaptive-epsilon | --fixed-epsilon\n"
        << "  --recursive-bisection --bisection-leaf-bins N\n"
        << "  --bisection-max-depth N\n"
        << "  --bisection-degree-limit N --bisection-fm-passes N\n"
        << "  --bisection-balance-tolerance X\n"
        << "  --bisection-position-seeded\n"
        << "  --bisection-surplus-only\n"
        << "  --bisection-leaf-curve-order\n"
        << "  --bisection-leaf-nearest-capacity\n"
        << "  --bisection-leaf-hpwl-guided\n"
        << "  --bisection-leaf-local-radius N\n"
        << "  --bisection-atomic-coarsening\n"
        << "  --bisection-atomic-max-nodes N --bisection-atomic-degree N\n"
        << "  --bisection-atomic-rounds N --bisection-atomic-release-depth N\n"
        << "  --bisection-axis-rank-map [--bisection-axis-map-only]\n"
        << "  --bisection-curve-rank-map [--bisection-curve-map-only]\n"
        << "  --coarse-capacity-flow --coarse-flow-bins N\n"
        << "  --coarse-flow-bins-x N --coarse-flow-bins-y N --coarse-flow-passes N\n"
        << "  --coarse-flow-exact-anchors\n"
        << "  --coarse-flow-max-source-bins N --coarse-flow-max-nodes-per-source N\n"
        << "  --coarse-flow-anchor-grid-stride N\n"
        << "  --coarse-flow-connected-commodities\n"
        << "  --coarse-flow-commodity-degree N --coarse-flow-commodity-max-nodes N\n"
        << "  --transport-rounds N --transport-max-moves N\n"
        << "  --transport-source-bins N --transport-lookahead N\n"
        << "  --transport-group-size N --transport-group-degree-limit N\n"
        << "  --transport-group-destination-radius N\n"
        << "  --transport-group-collective-target\n"
        << "  --transport-weighted-groups\n"
        << "  --transport-preserve-bin-offset\n"
        << "  --transport-rigid-group-displacement\n"
        << "  --transport-rigid-rounds N      (-1=all transport rounds)\n"
        << "  --transport-atomic-groups\n"
        << "  --transport-atomic-source-active\n"
        << "  --transport-atomic-rounds N     (-1=all transport rounds)\n"
        << "  --atomic-componentwise-capacity\n"
        << "  --atomic-capacity-candidates N\n"
        << "  --atomic-source-batch\n"
        << "  --atomic-min-gain-ratio X --atomic-fallback-individual\n"
        << "  --identity-exchange-passes N --identity-exchange-candidates N\n"
        << "  --identity-exchange-radius N\n"
        << "  --transport-destination nearest|auction|hilbert\n"
        << "  --transport-global-hpwl-budget X\n"
        << "  --auction-candidates N --auction-shortlist N\n"
        << "  --auction-hpwl-degree-limit N --auction-density-price-weight X\n"
        << "  --auction-rounds N             (-1=all transport rounds)\n"
        << "  --transport-connectivity-order --transport-group-order\n"
        << "  --connectivity-degree-limit N\n"
        << "  --hilbert-rounds N --hilbert-window N\n"
        << "  --recovery-sweeps N --recovery-step-bins X\n"
        << "  --recovery-line-search N --recovery-degree-limit N\n"
        << "  --recovery-axis-separated\n"
        << "  --recovery-breakpoint-oracle\n"
        << "  --recovery-boundary-active [--recovery-boundary-epsilon-bins X]\n"
        << "  --recovery-dead-zone-escape\n"
        << "  --recovery-controlled-exact-noise [--recovery-noise-seed N]\n"
        << "  --recovery-compact-directions\n"
        << "  --recovery-net-batch-weight X --recovery-net-batch-degree-limit N\n"
        << "  --recovery-net-blocks [--recovery-net-block-only]\n"
        << "  --recovery-net-block-density-direction\n"
        << "  --recovery-net-block-contraction\n"
        << "  --recovery-net-block-degree-limit N\n"
        << "  --recovery-net-block-max-nodes N --recovery-net-block-max-blocks N\n"
        << "  --recovery-overflow-cap X      (-1=each move non-increasing)\n"
        << "  --recovery-objective-density-weight X\n"
        << "  --density-coordinate-sweeps N --density-coordinate-step-bins X\n"
        << "  --density-coordinate-line-search N\n"
        << "  --density-coordinate-no-node-moves\n"
        << "  --density-coordinate-hpwl-budget X\n"
        << "  --density-coordinate-overflow-band X\n"
        << "  --density-coordinate-net-blocks\n"
        << "  --density-coordinate-capacity-assignment\n"
        << "  --density-coordinate-assignment-bins N\n"
        << "  --density-coordinate-hierarchical-clusters\n"
        << "  --density-coordinate-elastic-clusters\n"
        << "  --density-coordinate-cooperative-batch\n"
        << "  --density-coordinate-heavy-edge-growth\n"
        << "  --density-coordinate-adjacent-bin-only\n"
        << "  --density-coordinate-adjacent-bin-radius N\n"
        << "  --density-coordinate-hierarchy-levels N\n"
        << "  --density-coordinate-cluster-max-nodes N\n"
        << "  --density-coordinate-cluster-candidates N\n"
        << "  --density-coordinate-cluster-max-clusters N\n"
        << "  --density-coordinate-block-degree N\n"
        << "  --density-coordinate-block-max-nodes N\n"
        << "  --density-coordinate-block-max-blocks N\n"
        << "  --compact-sweeps N --compact-step-bins X\n"
        << "  --compact-line-search N --compact-overflow-cap X\n"
        << "  --swap-sweeps N --swap-radius-bins N --swap-candidates N\n"
        << "  --swap-exact-shortlist N --swap-degree-limit N --swap-net-aware\n"
        << "  --swap-density-guided --swap-density-hpwl-budget X\n"
        << "  --swap-density-score-weight X\n"
        << "  --swap-allow-unequal-shapes\n"
        << "  --shape-permutation-sweeps N\n"
        << "  --post-recovery-sweeps N --post-swap-sweeps N\n"
        << "  --anchor-assignment-sweeps N --assignment-radius-bins N\n"
        << "  --assignment-candidates N --assignment-max-bids-per-node N\n"
        << "  --assignment-epsilon X\n"
        << "  --snapshot-every N --seed N --sigma-ratio X --log-every N\n";
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_help();
            std::exit(0);
        } else if (arg == "--benchmark") {
            options.benchmark = value_after(i, argc, argv);
        } else if (arg == "--output") {
            options.output = value_after(i, argc, argv);
        } else if (arg == "--initial-placement") {
            options.initial_placement = value_after(i, argc, argv);
        } else if (arg == "--bins") {
            options.place.bins_x = options.place.bins_y =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bins-x") {
            options.place.bins_x = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bins-y") {
            options.place.bins_y = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--iterations") {
            options.place.iterations = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--multilevel") {
            options.multilevel = true;
        } else if (arg == "--level-iterations") {
            options.level_iterations = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--threads") {
            options.place.threads = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--optimizer") {
            options.place.optimizer = ea::parse_optimizer(value_after(i, argc, argv));
        } else if (arg == "--lambda-policy") {
            options.place.lambda.policy = ea::parse_lambda_policy(value_after(i, argc, argv));
        } else if (arg == "--hpwl-epsilon") {
            options.place.hpwl_epsilon = std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-epsilon") {
            options.place.density_epsilon = std::stod(value_after(i, argc, argv));
        } else if (arg == "--active-power") {
            options.place.active_power = std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-active-power") {
            options.place.density_active_power = std::stod(value_after(i, argc, argv));
        } else if (arg == "--step-fraction") {
            options.place.step_fraction = std::stod(value_after(i, argc, argv));
        } else if (arg == "--late-switch-iteration") {
            options.place.late_stage.switch_iteration =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--late-step-multiplier") {
            options.place.late_stage.step_multiplier =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--late-lambda-multiplier") {
            options.place.late_stage.lambda_multiplier =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--late-reset-optimizer") {
            options.place.late_stage.reset_optimizer = true;
        } else if (arg == "--max-step-multiplier") {
            options.place.max_step_multiplier = std::stod(value_after(i, argc, argv));
        } else if (arg == "--batch-backtracking") {
            options.place.batch_acceptance.enabled = true;
        } else if (arg == "--batch-infeasible") {
            options.place.batch_acceptance.allow_infeasible = true;
        } else if (arg == "--batch-hpwl-budget") {
            options.place.batch_acceptance.max_infeasible_hpwl_increase =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--batch-trials") {
            options.place.batch_acceptance.max_trials =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--batch-shrink") {
            options.place.batch_acceptance.shrink =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--beta1") {
            options.place.beta1 = std::stod(value_after(i, argc, argv));
        } else if (arg == "--beta2") {
            options.place.beta2 = std::stod(value_after(i, argc, argv));
        } else if (arg == "--momentum") {
            options.place.momentum = std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-weight-scale") {
            options.place.lambda.density_weight_scale =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--net-batch-weight") {
            options.place.net_batch.weight =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--net-batch-degree-limit") {
            options.place.net_batch.degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-net-share") {
            options.place.density_net_share.weight =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-net-degree-limit") {
            options.place.density_net_share.degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-net-zero-threshold") {
            options.place.density_net_share.zero_gradient_threshold =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-hpwl-orthogonal-projection") {
            options.place.density_net_share.hpwl_orthogonal_projection =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-net-share-until") {
            options.place.density_net_share.until_iteration =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--regional-price") {
            options.place.regional_price.enabled = true;
        } else if (arg == "--regional-price-interval") {
            options.place.regional_price.update_interval =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--regional-price-rho") {
            options.place.regional_price.rho =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--regional-price-target") {
            options.place.regional_price.target_excess =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--regional-price-min") {
            options.place.regional_price.min_price =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--regional-price-max") {
            options.place.regional_price.max_price =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--reference-hpwl-delta") {
            options.place.lambda.reference_hpwl_delta =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--lambda-update-interval") {
            options.place.lambda.update_interval =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--lambda-trajectory-horizon") {
            options.place.lambda.trajectory_horizon =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--target-density") {
            options.place.target_density = std::stod(value_after(i, argc, argv));
        } else if (arg == "--stop-overflow") {
            options.place.stop_overflow = std::stod(value_after(i, argc, argv));
        } else if (arg == "--degree-limit") {
            options.place.degree_limit = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--adaptive-epsilon") {
            options.place.adaptive_epsilon = true;
        } else if (arg == "--fixed-epsilon") {
            options.place.adaptive_epsilon = false;
        } else if (arg == "--adaptive-interval") {
            options.place.adaptive_interval = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--adaptive-window") {
            options.place.adaptive_window = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--adaptive-gain") {
            options.place.adaptive_gain = std::stod(value_after(i, argc, argv));
        } else if (arg == "--adaptive-deadband") {
            options.place.adaptive_deadband = std::stod(value_after(i, argc, argv));
        } else if (arg == "--adaptive-min-scale") {
            options.place.adaptive_min_scale = std::stod(value_after(i, argc, argv));
        } else if (arg == "--adaptive-max-scale") {
            options.place.adaptive_max_scale = std::stod(value_after(i, argc, argv));
        } else if (arg == "--adaptive-max-log-step") {
            options.place.adaptive_max_log_step = std::stod(value_after(i, argc, argv));
        } else if (arg == "--lower-overflow") {
            options.place.lower_overflow = std::stod(value_after(i, argc, argv));
        } else if (arg == "--recursive-bisection") {
            options.place.bisection.enabled = true;
        } else if (arg == "--bisection-leaf-bins") {
            options.place.bisection.leaf_bins =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bisection-max-depth") {
            options.place.bisection.max_depth =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bisection-degree-limit") {
            options.place.bisection.degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bisection-fm-passes") {
            options.place.bisection.fm_passes =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bisection-balance-tolerance") {
            options.place.bisection.balance_tolerance =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--bisection-position-seeded") {
            options.place.bisection.position_seeded = true;
        } else if (arg == "--bisection-surplus-only") {
            options.place.bisection.surplus_only = true;
        } else if (arg == "--bisection-leaf-curve-order") {
            options.place.bisection.leaf_curve_order = true;
        } else if (arg == "--bisection-leaf-nearest-capacity") {
            options.place.bisection.leaf_nearest_capacity = true;
        } else if (arg == "--bisection-leaf-hpwl-guided") {
            options.place.bisection.leaf_hpwl_guided = true;
        } else if (arg == "--bisection-leaf-local-radius") {
            options.place.bisection.leaf_local_radius =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bisection-atomic-coarsening") {
            options.place.bisection.atomic_coarsening = true;
        } else if (arg == "--bisection-atomic-max-nodes") {
            options.place.bisection.atomic_max_nodes =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bisection-atomic-degree") {
            options.place.bisection.atomic_degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bisection-atomic-rounds") {
            options.place.bisection.atomic_rounds =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bisection-atomic-release-depth") {
            options.place.bisection.atomic_release_depth =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--bisection-axis-rank-map") {
            options.place.bisection.axis_rank_map = true;
        } else if (arg == "--bisection-axis-map-only") {
            options.place.bisection.axis_rank_map = true;
            options.place.bisection.axis_map_only = true;
        } else if (arg == "--bisection-curve-rank-map") {
            options.place.bisection.curve_rank_map = true;
        } else if (arg == "--bisection-curve-map-only") {
            options.place.bisection.curve_rank_map = true;
            options.place.bisection.curve_map_only = true;
        } else if (arg == "--coarse-capacity-flow") {
            options.place.coarse_flow.enabled = true;
        } else if (arg == "--coarse-flow-bins") {
            options.place.coarse_flow.bins_x = options.place.coarse_flow.bins_y =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-bins-x") {
            options.place.coarse_flow.bins_x =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-bins-y") {
            options.place.coarse_flow.bins_y =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-passes") {
            options.place.coarse_flow.passes =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-exact-anchors") {
            options.place.coarse_flow.exact_anchors = true;
        } else if (arg == "--coarse-flow-max-source-bins") {
            options.place.coarse_flow.max_source_bins =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-max-nodes-per-source") {
            options.place.coarse_flow.max_nodes_per_source =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-anchor-grid-stride") {
            options.place.coarse_flow.anchor_grid_stride =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-connected-commodities") {
            options.place.coarse_flow.connected_commodities = true;
        } else if (arg == "--coarse-flow-commodity-multiscale") {
            options.place.coarse_flow.connected_commodities = true;
            options.place.coarse_flow.commodity_multiscale = true;
        } else if (arg == "--coarse-flow-commodity-axis-split") {
            options.place.coarse_flow.connected_commodities = true;
            options.place.coarse_flow.commodity_multiscale = true;
            options.place.coarse_flow.commodity_axis_split = true;
        } else if (arg == "--coarse-flow-commodity-hpwl-guard") {
            options.place.coarse_flow.connected_commodities = true;
            options.place.coarse_flow.commodity_hpwl_guard = true;
        } else if (arg == "--coarse-flow-joint-commodities") {
            options.place.coarse_flow.connected_commodities = true;
            options.place.coarse_flow.commodity_multiscale = true;
            options.place.coarse_flow.commodity_axis_split = true;
            options.place.coarse_flow.joint_assignment = true;
        } else if (arg == "--coarse-flow-joint-deferred") {
            options.place.coarse_flow.connected_commodities = true;
            options.place.coarse_flow.commodity_multiscale = true;
            options.place.coarse_flow.commodity_axis_split = true;
            options.place.coarse_flow.joint_assignment = true;
            options.place.coarse_flow.joint_allow_nonmonotone = true;
        } else if (arg == "--coarse-flow-joint-proposals") {
            options.place.coarse_flow.joint_proposals_per_edge =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-joint-no-edge-area-filter") {
            options.place.coarse_flow.connected_commodities = true;
            options.place.coarse_flow.commodity_multiscale = true;
            options.place.coarse_flow.commodity_axis_split = true;
            options.place.coarse_flow.joint_assignment = true;
            options.place.coarse_flow.joint_allow_nonmonotone = true;
            options.place.coarse_flow.joint_ignore_edge_area = true;
        } else if (arg == "--coarse-flow-joint-hpwl-budget") {
            options.place.coarse_flow.joint_hpwl_budget_ratio =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-commodity-degree") {
            options.place.coarse_flow.commodity_degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--coarse-flow-commodity-max-nodes") {
            options.place.coarse_flow.commodity_max_nodes =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-rounds") {
            options.place.transport.rounds = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-max-moves") {
            options.place.transport.max_moves = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-source-bins") {
            options.place.transport.max_source_bins = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-lookahead") {
            options.place.transport.candidate_lookahead = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-group-size") {
            options.place.transport.group_size = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-group-degree-limit") {
            options.place.transport.group_degree_limit = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-group-destination-radius") {
            options.place.transport.group_destination_radius =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-group-collective-target") {
            options.place.transport.group_collective_target = true;
        } else if (arg == "--transport-weighted-groups") {
            options.place.transport.weighted_grouping = true;
        } else if (arg == "--transport-preserve-bin-offset") {
            options.place.transport.preserve_bin_offset = true;
        } else if (arg == "--transport-rigid-group-displacement") {
            options.place.transport.rigid_group_displacement = true;
        } else if (arg == "--transport-rigid-rounds") {
            options.place.transport.rigid_rounds = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-atomic-groups") {
            options.place.transport.atomic_groups = true;
        } else if (arg == "--transport-atomic-source-active") {
            options.place.transport.atomic_groups = true;
            options.place.transport.atomic_source_active = true;
        } else if (arg == "--transport-atomic-rounds") {
            options.place.transport.atomic_rounds = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--atomic-componentwise-capacity") {
            options.place.transport.atomic_componentwise_capacity = true;
        } else if (arg == "--atomic-capacity-candidates") {
            options.place.transport.atomic_capacity_candidates =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--atomic-source-batch") {
            options.place.transport.atomic_source_batch = true;
        } else if (arg == "--atomic-min-gain-ratio") {
            options.place.transport.atomic_min_gain_ratio =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--atomic-fallback-individual") {
            options.place.transport.atomic_fallback_individual = true;
        } else if (arg == "--identity-exchange-passes") {
            options.place.transport.identity_exchange_passes =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--identity-exchange-candidates") {
            options.place.transport.identity_exchange_candidates =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--identity-exchange-radius") {
            options.place.transport.identity_exchange_radius =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-destination") {
            options.place.transport.destination_mode =
                ea::parse_transport_destination_mode(value_after(i, argc, argv));
        } else if (arg == "--transport-global-hpwl-budget") {
            options.place.transport.global_hpwl_budget =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--auction-candidates") {
            options.place.transport.auction_candidates = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--auction-shortlist") {
            options.place.transport.auction_shortlist = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--auction-hpwl-degree-limit") {
            options.place.transport.auction_hpwl_degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--auction-density-price-weight") {
            options.place.transport.auction_density_price_weight =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--auction-rounds") {
            options.place.transport.auction_rounds = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--transport-connectivity-order") {
            options.place.transport.connectivity_order = true;
        } else if (arg == "--transport-group-order") {
            options.place.transport.connectivity_order = false;
        } else if (arg == "--connectivity-degree-limit") {
            options.place.transport.connectivity_degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--hilbert-rounds") {
            options.place.transport.hilbert_rounds = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--hilbert-window") {
            options.place.transport.hilbert_window = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--recovery-sweeps") {
            options.place.recovery.sweeps = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--recovery-step-bins") {
            options.place.recovery.step_bins = std::stod(value_after(i, argc, argv));
        } else if (arg == "--recovery-line-search") {
            options.place.recovery.line_search_steps = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--recovery-degree-limit") {
            options.place.recovery.degree_limit = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--recovery-axis-separated") {
            options.place.recovery.axis_separated = true;
        } else if (arg == "--recovery-breakpoint-oracle") {
            options.place.recovery.breakpoint_oracle = true;
        } else if (arg == "--recovery-boundary-active") {
            options.place.recovery.boundary_active = true;
        } else if (arg == "--recovery-boundary-epsilon-bins") {
            options.place.recovery.boundary_epsilon_bins =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--recovery-dead-zone-escape") {
            options.place.recovery.dead_zone_escape = true;
        } else if (arg == "--recovery-controlled-exact-noise") {
            options.place.recovery.controlled_exact_noise = true;
        } else if (arg == "--recovery-noise-seed") {
            options.place.recovery.exact_noise_seed =
                std::stoull(value_after(i, argc, argv));
        } else if (arg == "--recovery-compact-directions") {
            options.place.recovery.compact_directions = true;
        } else if (arg == "--recovery-net-batch-weight") {
            options.place.recovery.net_batch_weight =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--recovery-net-batch-degree-limit") {
            options.place.recovery.net_batch_degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--recovery-net-blocks") {
            options.place.recovery.net_block = true;
        } else if (arg == "--recovery-net-block-only") {
            options.place.recovery.net_block = true;
            options.place.recovery.node_moves = false;
        } else if (arg == "--recovery-net-block-density-direction") {
            options.place.recovery.net_block = true;
            options.place.recovery.net_block_density_direction = true;
        } else if (arg == "--recovery-net-block-contraction") {
            options.place.recovery.net_block = true;
            options.place.recovery.net_block_contraction = true;
        } else if (arg == "--recovery-net-block-degree-limit") {
            options.place.recovery.net_block_degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--recovery-net-block-max-nodes") {
            options.place.recovery.net_block_max_nodes =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--recovery-net-block-max-blocks") {
            options.place.recovery.net_block_max_blocks =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--recovery-hpwl-epsilon") {
            options.place.recovery.hpwl_epsilon = std::stod(value_after(i, argc, argv));
        } else if (arg == "--recovery-overflow-cap") {
            options.place.recovery.overflow_cap = std::stod(value_after(i, argc, argv));
        } else if (arg == "--recovery-objective-density-weight") {
            options.place.recovery.objective_density_weight =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-sweeps") {
            options.place.density_coordinate.sweeps =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-step-bins") {
            options.place.density_coordinate.step_bins =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-line-search") {
            options.place.density_coordinate.line_search_steps =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-no-node-moves") {
            options.place.density_coordinate.node_moves = false;
        } else if (arg == "--density-coordinate-hpwl-budget") {
            options.place.density_coordinate.hpwl_budget_fraction =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-overflow-band") {
            options.place.density_coordinate.overflow_band_fraction =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-net-blocks") {
            options.place.density_coordinate.net_blocks = true;
        } else if (arg == "--density-coordinate-capacity-assignment") {
            options.place.density_coordinate.capacity_assignment = true;
            options.place.density_coordinate.net_blocks = true;
        } else if (arg == "--density-coordinate-assignment-bins") {
            options.place.density_coordinate.assignment_bins =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-hierarchical-clusters") {
            options.place.density_coordinate.hierarchical_clusters = true;
        } else if (arg == "--density-coordinate-elastic-clusters") {
            options.place.density_coordinate.elastic_clusters = true;
            options.place.density_coordinate.capacity_assignment = true;
            options.place.density_coordinate.net_blocks = true;
        } else if (arg == "--density-coordinate-cooperative-batch") {
            options.place.density_coordinate.cooperative_batch = true;
            options.place.density_coordinate.hierarchical_clusters = true;
            options.place.density_coordinate.capacity_assignment = true;
            options.place.density_coordinate.net_blocks = true;
        } else if (arg == "--density-coordinate-heavy-edge-growth") {
            options.place.density_coordinate.heavy_edge_growth = true;
        } else if (arg == "--density-coordinate-adjacent-bin-only") {
            options.place.density_coordinate.adjacent_bin_only = true;
        } else if (arg == "--density-coordinate-adjacent-bin-radius") {
            options.place.density_coordinate.adjacent_bin_radius =
                std::stoi(value_after(i, argc, argv));
            options.place.density_coordinate.hierarchical_clusters = true;
            options.place.density_coordinate.capacity_assignment = true;
            options.place.density_coordinate.net_blocks = true;
        } else if (arg == "--density-coordinate-hierarchy-levels") {
            options.place.density_coordinate.hierarchy_levels =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-cluster-max-nodes") {
            options.place.density_coordinate.cluster_max_nodes =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-cluster-candidates") {
            options.place.density_coordinate.cluster_candidate_bins =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-cluster-max-clusters") {
            options.place.density_coordinate.cluster_max_clusters =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-block-degree") {
            options.place.density_coordinate.block_degree_limit =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-block-max-nodes") {
            options.place.density_coordinate.block_max_nodes =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--density-coordinate-block-max-blocks") {
            options.place.density_coordinate.block_max_blocks =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--compact-sweeps") {
            options.place.compact_recovery.sweeps =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--compact-step-bins") {
            options.place.compact_recovery.max_step_bins =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--compact-line-search") {
            options.place.compact_recovery.line_search_steps =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--compact-overflow-cap") {
            options.place.compact_recovery.overflow_cap =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--compact-snapshot-every") {
            options.place.compact_recovery.snapshot_every =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--swap-sweeps") {
            options.place.swap_recovery.sweeps = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--swap-radius-bins") {
            options.place.swap_recovery.radius_bins = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--swap-candidates") {
            options.place.swap_recovery.candidates = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--swap-exact-shortlist") {
            options.place.swap_recovery.exact_shortlist =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--swap-net-aware") {
            options.place.swap_recovery.net_aware = true;
        } else if (arg == "--swap-density-guided") {
            options.place.swap_recovery.density_guided = true;
        } else if (arg == "--swap-density-hpwl-budget") {
            options.place.swap_recovery.density_hpwl_budget =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--swap-density-score-weight") {
            options.place.swap_recovery.density_score_weight =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--swap-allow-unequal-shapes") {
            options.place.swap_recovery.allow_unequal_shapes = true;
        } else if (arg == "--swap-degree-limit") {
            options.place.swap_recovery.degree_limit = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--shape-permutation-sweeps") {
            options.place.swap_recovery.permutation_sweeps =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--post-recovery-sweeps") {
            options.place.post_recovery.sweeps = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--post-swap-sweeps") {
            options.place.post_swap_recovery.sweeps = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--anchor-assignment-sweeps") {
            options.place.assignment_recovery.assignment_sweeps =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--assignment-radius-bins") {
            options.place.assignment_recovery.assignment_radius_bins =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--assignment-candidates") {
            options.place.assignment_recovery.assignment_candidates =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--assignment-max-bids-per-node") {
            options.place.assignment_recovery.assignment_max_bids_per_node =
                std::stoi(value_after(i, argc, argv));
        } else if (arg == "--assignment-epsilon") {
            options.place.assignment_recovery.assignment_epsilon =
                std::stod(value_after(i, argc, argv));
        } else if (arg == "--snapshot-every") {
            options.place.snapshot_every = std::stoi(value_after(i, argc, argv));
        } else if (arg == "--seed") {
            options.seed = std::stoull(value_after(i, argc, argv));
        } else if (arg == "--sigma-ratio") {
            options.sigma_ratio = std::stod(value_after(i, argc, argv));
        } else if (arg == "--log-every") {
            options.place.log_every = std::stoi(value_after(i, argc, argv));
        } else {
            throw std::invalid_argument("unknown option: " + arg);
        }
    }
    if (options.benchmark.empty() || options.output.empty()) {
        throw std::invalid_argument("--benchmark and --output are required");
    }
    options.place.output_dir = options.output;
    options.place.compact_recovery.output_dir = options.output;
    options.place.lambda.stop_overflow = options.place.stop_overflow;
    if (options.place.post_recovery.sweeps > 0) {
        const int sweeps = options.place.post_recovery.sweeps;
        options.place.post_recovery = options.place.recovery;
        options.place.post_recovery.sweeps = sweeps;
    }
    if (options.place.post_swap_recovery.sweeps > 0) {
        const int sweeps = options.place.post_swap_recovery.sweeps;
        options.place.post_swap_recovery = options.place.swap_recovery;
        options.place.post_swap_recovery.sweeps = sweeps;
        options.place.post_swap_recovery.permutation_sweeps = 0;
    }
    return options;
}

void write_summary(const Options& options, const ea::Database& db,
                   const ea::PlaceResult& result) {
    std::ofstream output(options.output / "summary.txt");
    if (!output) throw std::runtime_error("cannot create summary.txt");
    output << std::setprecision(12);
    output << "benchmark=" << db.benchmark_base << '\n';
    output << "initial_placement=" << options.initial_placement.string() << '\n';
    output << "initialization_mode="
           << (options.initial_placement.empty()
                   ? "fresh-center-gaussian"
                   : "explicit-stage-continuation")
           << '\n';
    output << "initialization_seed=" << options.seed << '\n';
    output << "initialization_sigma_ratio=" << options.sigma_ratio << '\n';
    output << "benchmark_raw_pl=" << db.raw_pl_path << '\n';
    output << "model=exact-hpwl-plus-exact-grid-overlap\n";
    output << "legalization=disabled\n";
    output << "bins=" << options.place.bins_x << 'x' << options.place.bins_y << '\n';
    output << "multilevel=" << (options.multilevel ? "true" : "false") << '\n';
    output << "level_iterations=" << options.level_iterations << '\n';
    output << "iterations=" << options.place.iterations << '\n';
    output << "threads_requested=" << options.place.threads << '\n';
    output << "threads_effective=" << omp_get_max_threads() << '\n';
    output << "optimizer=" << ea::optimizer_name(options.place.optimizer) << '\n';
    output << "lambda_policy=" << ea::lambda_policy_name(options.place.lambda.policy) << '\n';
    output << "density_weight_scale=" << options.place.lambda.density_weight_scale << '\n';
    output << "net_batch_weight=" << options.place.net_batch.weight << '\n';
    output << "net_batch_degree_limit=" << options.place.net_batch.degree_limit << '\n';
    output << "density_net_share_weight=" << options.place.density_net_share.weight << '\n';
    output << "density_net_share_degree_limit=" << options.place.density_net_share.degree_limit << '\n';
    output << "density_net_share_zero_threshold="
           << options.place.density_net_share.zero_gradient_threshold << '\n';
    output << "density_hpwl_orthogonal_projection="
           << options.place.density_net_share.hpwl_orthogonal_projection << '\n';
    output << "density_net_share_until="
           << options.place.density_net_share.until_iteration << '\n';
    output << "regional_price_enabled="
           << (options.place.regional_price.enabled ? "true" : "false") << '\n';
    output << "regional_price_interval="
           << options.place.regional_price.update_interval << '\n';
    output << "regional_price_rho=" << options.place.regional_price.rho << '\n';
    output << "step_fraction=" << options.place.step_fraction << '\n';
    output << "late_switch_iteration="
           << options.place.late_stage.switch_iteration << '\n';
    output << "late_step_multiplier="
           << options.place.late_stage.step_multiplier << '\n';
    output << "late_lambda_multiplier="
           << options.place.late_stage.lambda_multiplier << '\n';
    output << "late_reset_optimizer="
           << (options.place.late_stage.reset_optimizer ? "true" : "false") << '\n';
    output << "batch_backtracking="
           << (options.place.batch_acceptance.enabled ? "true" : "false") << '\n';
    output << "batch_allow_infeasible="
           << (options.place.batch_acceptance.allow_infeasible ? "true" : "false") << '\n';
    output << "batch_infeasible_hpwl_budget="
           << options.place.batch_acceptance.max_infeasible_hpwl_increase << '\n';
    output << "batch_max_trials=" << options.place.batch_acceptance.max_trials << '\n';
    output << "batch_shrink=" << options.place.batch_acceptance.shrink << '\n';
    output << "accepted_batches=" << result.accepted_batches << '\n';
    output << "rejected_batches=" << result.rejected_batches << '\n';
    output << "batch_trials=" << result.batch_trials << '\n';
    output << "active_power=" << options.place.active_power << '\n';
    output << "density_active_power=" << options.place.density_active_power << '\n';
    output << "target_density=" << options.place.target_density << '\n';
    output << "stop_overflow=" << options.place.stop_overflow << '\n';
    output << "hpwl_epsilon=" << options.place.hpwl_epsilon << '\n';
    output << "density_epsilon=" << options.place.density_epsilon << '\n';
    output << "adaptive_epsilon=" << (options.place.adaptive_epsilon ? "true" : "false") << '\n';
    output << "recursive_bisection="
           << (options.place.bisection.enabled ? "true" : "false") << '\n';
    output << "bisection_leaf_bins=" << options.place.bisection.leaf_bins << '\n';
    output << "bisection_max_depth=" << options.place.bisection.max_depth << '\n';
    output << "bisection_degree_limit="
           << options.place.bisection.degree_limit << '\n';
    output << "bisection_fm_passes=" << options.place.bisection.fm_passes << '\n';
    output << "bisection_balance_tolerance="
           << options.place.bisection.balance_tolerance << '\n';
    output << "bisection_position_seeded="
           << (options.place.bisection.position_seeded ? "true" : "false") << '\n';
    output << "bisection_surplus_only="
           << (options.place.bisection.surplus_only ? "true" : "false") << '\n';
    output << "bisection_leaf_curve_order="
           << (options.place.bisection.leaf_curve_order ? "true" : "false") << '\n';
    output << "bisection_leaf_nearest_capacity="
           << (options.place.bisection.leaf_nearest_capacity ? "true" : "false") << '\n';
    output << "bisection_leaf_hpwl_guided="
           << (options.place.bisection.leaf_hpwl_guided ? "true" : "false") << '\n';
    output << "bisection_leaf_local_radius="
           << options.place.bisection.leaf_local_radius << '\n';
    output << "bisection_atomic_coarsening="
           << (options.place.bisection.atomic_coarsening ? "true" : "false")
           << '\n';
    output << "bisection_atomic_max_nodes="
           << options.place.bisection.atomic_max_nodes << '\n';
    output << "bisection_atomic_degree="
           << options.place.bisection.atomic_degree_limit << '\n';
    output << "bisection_atomic_rounds="
           << options.place.bisection.atomic_rounds << '\n';
    output << "bisection_atomic_release_depth="
           << options.place.bisection.atomic_release_depth << '\n';
    output << "bisection_axis_rank_map="
           << (options.place.bisection.axis_rank_map ? "true" : "false") << '\n';
    output << "bisection_axis_map_only="
           << (options.place.bisection.axis_map_only ? "true" : "false") << '\n';
    output << "bisection_curve_rank_map="
           << (options.place.bisection.curve_rank_map ? "true" : "false") << '\n';
    output << "bisection_curve_map_only="
           << (options.place.bisection.curve_map_only ? "true" : "false") << '\n';
    output << "bisection_splits=" << result.bisection.splits << '\n';
    output << "bisection_leaves=" << result.bisection.leaves << '\n';
    output << "bisection_max_depth=" << result.bisection.max_depth << '\n';
    output << "bisection_fm_moves=" << result.bisection.fm_moves << '\n';
    output << "bisection_atomic_groups=" << result.bisection.atomic_groups << '\n';
    output << "bisection_atomic_merges=" << result.bisection.atomic_merges << '\n';
    output << "bisection_atomic_unit_moves="
           << result.bisection.atomic_unit_moves << '\n';
    output << "bisection_cut_net_weight=" << result.bisection.cut_net_weight << '\n';
    output << "bisection_capacity_error="
           << result.bisection.maximum_capacity_error << '\n';
    output << "bisection_final_hpwl=" << result.bisection.final_hpwl << '\n';
    output << "bisection_final_overflow=" << result.bisection.final_overflow << '\n';
    output << "bisection_seconds=" << result.bisection.wall_seconds << '\n';
    output << "coarse_capacity_flow="
           << (options.place.coarse_flow.enabled ? "true" : "false") << '\n';
    output << "coarse_flow_bins=" << options.place.coarse_flow.bins_x << 'x'
           << options.place.coarse_flow.bins_y << '\n';
    output << "coarse_flow_passes=" << result.coarse_flow.passes << '\n';
    output << "coarse_flow_exact_anchors="
           << (options.place.coarse_flow.exact_anchors ? "true" : "false") << '\n';
    output << "coarse_flow_max_source_bins="
           << options.place.coarse_flow.max_source_bins << '\n';
    output << "coarse_flow_max_nodes_per_source="
           << options.place.coarse_flow.max_nodes_per_source << '\n';
    output << "coarse_flow_anchor_grid_stride="
           << options.place.coarse_flow.anchor_grid_stride << '\n';
    output << "coarse_flow_connected_commodities="
           << (options.place.coarse_flow.connected_commodities ? "true" : "false") << '\n';
    output << "coarse_flow_commodity_degree="
           << options.place.coarse_flow.commodity_degree_limit << '\n';
    output << "coarse_flow_commodity_max_nodes="
           << options.place.coarse_flow.commodity_max_nodes << '\n';
    output << "coarse_flow_commodity_multiscale="
           << (options.place.coarse_flow.commodity_multiscale ? "true" : "false") << '\n';
    output << "coarse_flow_commodity_axis_split="
           << (options.place.coarse_flow.commodity_axis_split ? "true" : "false") << '\n';
    output << "coarse_flow_commodity_hpwl_guard="
           << (options.place.coarse_flow.commodity_hpwl_guard ? "true" : "false") << '\n';
    output << "coarse_flow_joint_assignment="
           << (options.place.coarse_flow.joint_assignment ? "true" : "false") << '\n';
    output << "coarse_flow_joint_allow_nonmonotone="
           << (options.place.coarse_flow.joint_allow_nonmonotone ? "true" : "false") << '\n';
    output << "coarse_flow_joint_proposals_per_edge="
           << options.place.coarse_flow.joint_proposals_per_edge << '\n';
    output << "coarse_flow_joint_ignore_edge_area="
           << (options.place.coarse_flow.joint_ignore_edge_area ? "true" : "false") << '\n';
    output << "coarse_flow_joint_hpwl_budget_ratio="
           << options.place.coarse_flow.joint_hpwl_budget_ratio << '\n';
    output << "coarse_flow_edges=" << result.coarse_flow.flow_edges << '\n';
    output << "coarse_flow_moves=" << result.coarse_flow.moves << '\n';
    output << "coarse_flow_exact_anchor_trials="
           << result.coarse_flow.exact_anchor_trials << '\n';
    output << "coarse_flow_exact_anchor_rejections="
           << result.coarse_flow.exact_anchor_rejections << '\n';
    output << "coarse_flow_planned_area=" << result.coarse_flow.planned_area << '\n';
    output << "coarse_flow_moved_area=" << result.coarse_flow.moved_area << '\n';
    output << "coarse_flow_initial_coarse_overflow="
           << result.coarse_flow.initial_coarse_overflow << '\n';
    output << "coarse_flow_final_coarse_overflow="
           << result.coarse_flow.final_coarse_overflow << '\n';
    output << "coarse_flow_initial_hpwl=" << result.coarse_flow.initial_hpwl << '\n';
    output << "coarse_flow_initial_overflow="
           << result.coarse_flow.initial_overflow << '\n';
    output << "coarse_flow_final_hpwl=" << result.coarse_flow.final_hpwl << '\n';
    output << "coarse_flow_final_overflow="
           << result.coarse_flow.final_overflow << '\n';
    output << "coarse_flow_seconds=" << result.coarse_flow.wall_seconds << '\n';
    output << "coarse_flow_commodity_groups=" << result.coarse_flow.commodity_groups << '\n';
    output << "coarse_flow_commodity_moves=" << result.coarse_flow.commodity_moves << '\n';
    output << "coarse_flow_commodity_scale_trials="
           << result.coarse_flow.commodity_scale_trials << '\n';
    output << "coarse_flow_joint_batches=" << result.coarse_flow.joint_batches << '\n';
    output << "coarse_flow_joint_batch_trials="
           << result.coarse_flow.joint_batch_trials << '\n';
    output << "transport_rounds=" << result.transport.rounds << '\n';
    output << "transport_moves=" << result.transport.moves << '\n';
    output << "transport_groups=" << result.transport.groups << '\n';
    output << "transport_group_destination_radius="
           << options.place.transport.group_destination_radius << '\n';
    output << "transport_group_collective_target="
           << (options.place.transport.group_collective_target ? "true" : "false") << '\n';
    output << "transport_weighted_groups="
           << (options.place.transport.weighted_grouping ? "true" : "false") << '\n';
    output << "transport_preserve_bin_offset="
           << (options.place.transport.preserve_bin_offset ? "true" : "false") << '\n';
    output << "transport_rigid_group_displacement="
           << (options.place.transport.rigid_group_displacement ? "true" : "false") << '\n';
    output << "transport_rigid_rounds=" << options.place.transport.rigid_rounds << '\n';
    output << "transport_atomic_groups="
           << (options.place.transport.atomic_groups ? "true" : "false") << '\n';
    output << "transport_atomic_source_active="
           << (options.place.transport.atomic_source_active ? "true" : "false") << '\n';
    output << "transport_atomic_rounds=" << options.place.transport.atomic_rounds << '\n';
    output << "atomic_componentwise_capacity="
           << (options.place.transport.atomic_componentwise_capacity ? "true" : "false")
           << '\n';
    output << "atomic_capacity_candidates="
           << options.place.transport.atomic_capacity_candidates << '\n';
    output << "atomic_source_batch="
           << (options.place.transport.atomic_source_batch ? "true" : "false") << '\n';
    output << "atomic_min_gain_ratio="
           << options.place.transport.atomic_min_gain_ratio << '\n';
    output << "atomic_fallback_individual="
           << (options.place.transport.atomic_fallback_individual ? "true" : "false")
           << '\n';
    output << "identity_exchange_passes="
           << options.place.transport.identity_exchange_passes << '\n';
    output << "identity_exchange_candidates="
           << options.place.transport.identity_exchange_candidates << '\n';
    output << "identity_exchange_radius="
           << options.place.transport.identity_exchange_radius << '\n';
    output << "identity_exchange_attempts="
           << result.transport.identity_exchange_attempts << '\n';
    output << "identity_exchange_accepted="
           << result.transport.identity_exchange_accepted << '\n';
    output << "identity_exchange_permuted_nodes="
           << result.transport.identity_exchange_permuted_nodes << '\n';
    output << "identity_exchange_initial_hpwl="
           << result.transport.identity_initial_hpwl << '\n';
    output << "identity_exchange_final_hpwl="
           << result.transport.identity_final_hpwl << '\n';
    output << "identity_exchange_initial_overflow="
           << result.transport.identity_initial_overflow << '\n';
    output << "identity_exchange_final_overflow="
           << result.transport.identity_final_overflow << '\n';
    output << "identity_exchange_seconds="
           << result.transport.identity_wall_seconds << '\n';
    output << "transport_atomic_attempts=" << result.transport.atomic_attempts << '\n';
    output << "transport_atomic_accepted=" << result.transport.atomic_groups << '\n';
    output << "transport_atomic_fallbacks=" << result.transport.atomic_fallbacks << '\n';
    output << "transport_atomic_capacity_rejections="
           << result.transport.atomic_capacity_rejections << '\n';
    output << "transport_atomic_bid_trials=" << result.transport.atomic_bid_trials << '\n';
    output << "transport_atomic_feasible_bids="
           << result.transport.atomic_feasible_bids << '\n';
    output << "transport_atomic_rescued_groups="
           << result.transport.atomic_rescued_groups << '\n';
    output << "transport_atomic_batch_sources="
           << result.transport.atomic_batch_sources << '\n';
    output << "transport_atomic_batch_groups="
           << result.transport.atomic_batch_groups << '\n';
    output << "transport_atomic_batch_bids=" << result.transport.atomic_batch_bids << '\n';
    output << "transport_atomic_batch_rejections="
           << result.transport.atomic_batch_rejections << '\n';
    output << "transport_atomic_batch_commits="
           << result.transport.atomic_batch_commits << '\n';
    output << "transport_destination="
           << ea::transport_destination_mode_name(options.place.transport.destination_mode)
           << '\n';
    output << "transport_global_hpwl_budget="
           << options.place.transport.global_hpwl_budget << '\n';
    output << "auction_rounds=" << options.place.transport.auction_rounds << '\n';
    output << "transport_connectivity_order="
           << (options.place.transport.connectivity_order ? "true" : "false") << '\n';
    output << "hilbert_rounds=" << options.place.transport.hilbert_rounds << '\n';
    output << "hilbert_window=" << options.place.transport.hilbert_window << '\n';
    output << "transport_initial_hpwl=" << result.transport.initial_hpwl << '\n';
    output << "transport_initial_overflow=" << result.transport.initial_overflow << '\n';
    output << "transport_final_hpwl=" << result.transport.final_hpwl << '\n';
    output << "transport_final_overflow=" << result.transport.final_overflow << '\n';
    output << "transport_seconds=" << result.transport.wall_seconds << '\n';
    output << "density_coordinate_sweeps="
           << result.density_coordinate.sweeps << '\n';
    output << "density_coordinate_moves="
           << result.density_coordinate.moves << '\n';
    output << "density_coordinate_candidates="
           << result.density_coordinate.candidates << '\n';
    output << "density_coordinate_cluster_proposals="
           << result.density_coordinate.cluster_proposals << '\n';
    output << "density_coordinate_cluster_moves="
           << result.density_coordinate.cluster_moves << '\n';
    output << "density_coordinate_cooperative_batch_trials="
           << result.density_coordinate.cooperative_batch_trials << '\n';
    output << "density_coordinate_cooperative_batches="
           << result.density_coordinate.cooperative_batches << '\n';
    output << "density_coordinate_initial_hpwl="
           << result.density_coordinate.initial_hpwl << '\n';
    output << "density_coordinate_initial_overflow="
           << result.density_coordinate.initial_overflow << '\n';
    output << "density_coordinate_final_hpwl="
           << result.density_coordinate.final_hpwl << '\n';
    output << "density_coordinate_final_overflow="
           << result.density_coordinate.final_overflow << '\n';
    output << "density_coordinate_seconds="
           << result.density_coordinate.wall_seconds << '\n';
    output << "density_coordinate_hpwl_budget="
           << options.place.density_coordinate.hpwl_budget_fraction << '\n';
    output << "density_coordinate_overflow_band="
           << options.place.density_coordinate.overflow_band_fraction << '\n';
    output << "density_coordinate_hierarchical_clusters="
           << (options.place.density_coordinate.hierarchical_clusters
                   ? "true" : "false") << '\n';
    output << "density_coordinate_cooperative_batch="
           << (options.place.density_coordinate.cooperative_batch
                   ? "true" : "false") << '\n';
    output << "density_coordinate_heavy_edge_growth="
           << (options.place.density_coordinate.heavy_edge_growth
               ? "true" : "false") << '\n';
    output << "density_coordinate_adjacent_bin_only="
           << (options.place.density_coordinate.adjacent_bin_only
               ? "true" : "false") << '\n';
    output << "density_coordinate_adjacent_bin_radius="
           << options.place.density_coordinate.adjacent_bin_radius << '\n';
    output << "density_coordinate_hierarchy_levels="
           << options.place.density_coordinate.hierarchy_levels << '\n';
    output << "density_coordinate_cluster_max_nodes="
           << options.place.density_coordinate.cluster_max_nodes << '\n';
    output << "density_coordinate_cluster_candidates="
           << options.place.density_coordinate.cluster_candidate_bins << '\n';
    output << "density_coordinate_cluster_max_clusters="
           << options.place.density_coordinate.cluster_max_clusters << '\n';
    output << "recovery_sweeps=" << result.recovery.sweeps << '\n';
    output << "recovery_moves=" << result.recovery.moves << '\n';
    output << "recovery_initial_hpwl=" << result.recovery.initial_hpwl << '\n';
    output << "recovery_initial_overflow=" << result.recovery.initial_overflow << '\n';
    output << "recovery_final_hpwl=" << result.recovery.final_hpwl << '\n';
    output << "recovery_final_overflow=" << result.recovery.final_overflow << '\n';
    output << "recovery_seconds=" << result.recovery.wall_seconds << '\n';
    output << "recovery_overflow_cap=" << options.place.recovery.overflow_cap << '\n';
    output << "recovery_objective_density_weight="
           << options.place.recovery.objective_density_weight << '\n';
    output << "recovery_axis_separated="
           << (options.place.recovery.axis_separated ? "true" : "false") << '\n';
    output << "recovery_breakpoint_oracle="
           << (options.place.recovery.breakpoint_oracle ? "true" : "false") << '\n';
    output << "recovery_boundary_active="
           << (options.place.recovery.boundary_active ? "true" : "false") << '\n';
    output << "recovery_boundary_epsilon_bins="
           << options.place.recovery.boundary_epsilon_bins << '\n';
    output << "recovery_dead_zone_escape="
           << (options.place.recovery.dead_zone_escape ? "true" : "false") << '\n';
    output << "recovery_controlled_exact_noise="
           << (options.place.recovery.controlled_exact_noise ? "true" : "false") << '\n';
    output << "recovery_exact_noise_seed="
           << options.place.recovery.exact_noise_seed << '\n';
    output << "recovery_compact_directions="
           << (options.place.recovery.compact_directions ? "true" : "false") << '\n';
    output << "recovery_compact_direction_attempts="
           << result.recovery.compact_direction_attempts << '\n';
    output << "recovery_compact_direction_moves="
           << result.recovery.compact_direction_moves << '\n';
    output << "recovery_net_batch_weight="
           << options.place.recovery.net_batch_weight << '\n';
    output << "recovery_net_block="
           << (options.place.recovery.net_block ? "true" : "false") << '\n';
    output << "recovery_node_moves="
           << (options.place.recovery.node_moves ? "true" : "false") << '\n';
    output << "recovery_net_block_density_direction="
           << (options.place.recovery.net_block_density_direction
                   ? "true" : "false") << '\n';
    output << "recovery_net_block_contraction="
           << (options.place.recovery.net_block_contraction ? "true" : "false")
           << '\n';
    output << "recovery_net_block_degree_limit="
           << options.place.recovery.net_block_degree_limit << '\n';
    output << "recovery_net_block_max_nodes="
           << options.place.recovery.net_block_max_nodes << '\n';
    output << "recovery_net_block_max_blocks="
           << options.place.recovery.net_block_max_blocks << '\n';
    output << "recovery_net_block_attempts="
           << result.recovery.net_block_attempts << '\n';
    output << "recovery_net_block_moves="
           << result.recovery.net_block_moves << '\n';
    output << "recovery_net_block_nodes="
           << result.recovery.net_block_nodes << '\n';
    output << "recovery_net_block_contraction_attempts="
           << result.recovery.net_block_contraction_attempts << '\n';
    output << "recovery_net_block_contraction_moves="
           << result.recovery.net_block_contraction_moves << '\n';
    output << "recovery_net_block_contraction_nodes="
           << result.recovery.net_block_contraction_nodes << '\n';
    output << "compact_sweeps=" << result.compact_recovery.sweeps << '\n';
    output << "compact_candidates=" << result.compact_recovery.candidates << '\n';
    output << "compact_accepted_candidates="
           << result.compact_recovery.accepted_candidates << '\n';
    output << "compact_step_bins="
           << options.place.compact_recovery.max_step_bins << '\n';
    output << "compact_overflow_cap="
           << options.place.compact_recovery.overflow_cap << '\n';
    output << "compact_initial_hpwl="
           << result.compact_recovery.initial_hpwl << '\n';
    output << "compact_initial_overflow="
           << result.compact_recovery.initial_overflow << '\n';
    output << "compact_final_hpwl="
           << result.compact_recovery.final_hpwl << '\n';
    output << "compact_final_overflow="
           << result.compact_recovery.final_overflow << '\n';
    output << "compact_initial_spread="
           << result.compact_recovery.initial_compactness << '\n';
    output << "compact_final_spread="
           << result.compact_recovery.final_compactness << '\n';
    output << "compact_seconds="
           << result.compact_recovery.wall_seconds << '\n';
    output << "swap_recovery_sweeps=" << result.swap_recovery.sweeps << '\n';
    output << "swap_recovery_swaps=" << result.swap_recovery.swaps << '\n';
    output << "swap_recovery_initial_hpwl=" << result.swap_recovery.initial_hpwl << '\n';
    output << "swap_recovery_final_hpwl=" << result.swap_recovery.final_hpwl << '\n';
    output << "swap_recovery_initial_overflow=" << result.swap_recovery.initial_overflow << '\n';
    output << "swap_recovery_final_overflow=" << result.swap_recovery.final_overflow << '\n';
    output << "swap_recovery_seconds=" << result.swap_recovery.wall_seconds << '\n';
    output << "swap_exact_shortlist="
           << options.place.swap_recovery.exact_shortlist << '\n';
    output << "swap_density_guided="
           << (options.place.swap_recovery.density_guided ? "true" : "false") << '\n';
    output << "swap_density_hpwl_budget="
           << options.place.swap_recovery.density_hpwl_budget << '\n';
    output << "swap_density_score_weight="
           << options.place.swap_recovery.density_score_weight << '\n';
    output << "swap_allow_unequal_shapes="
           << (options.place.swap_recovery.allow_unequal_shapes ? "true" : "false") << '\n';
    output << "swap_net_aware="
           << (options.place.swap_recovery.net_aware ? "true" : "false") << '\n';
    output << "shape_permutation_sweeps="
           << result.swap_recovery.permutation_sweeps << '\n';
    output << "shape_permuted_nodes=" << result.swap_recovery.permuted_nodes << '\n';
    output << "post_recovery_sweeps=" << result.post_recovery.sweeps << '\n';
    output << "post_recovery_moves=" << result.post_recovery.moves << '\n';
    output << "post_recovery_final_hpwl=" << result.post_recovery.final_hpwl << '\n';
    output << "post_recovery_final_overflow=" << result.post_recovery.final_overflow << '\n';
    output << "post_recovery_seconds=" << result.post_recovery.wall_seconds << '\n';
    output << "post_swap_sweeps=" << result.post_swap_recovery.sweeps << '\n';
    output << "post_swap_swaps=" << result.post_swap_recovery.swaps << '\n';
    output << "post_swap_final_hpwl=" << result.post_swap_recovery.final_hpwl << '\n';
    output << "post_swap_final_overflow=" << result.post_swap_recovery.final_overflow << '\n';
    output << "post_swap_seconds=" << result.post_swap_recovery.wall_seconds << '\n';
    output << "assignment_sweeps="
           << result.assignment_recovery.assignment_sweeps << '\n';
    output << "assignment_radius_bins="
           << options.place.assignment_recovery.assignment_radius_bins << '\n';
    output << "assignment_candidates="
           << options.place.assignment_recovery.assignment_candidates << '\n';
    output << "assignment_max_bids_per_node="
           << options.place.assignment_recovery.assignment_max_bids_per_node << '\n';
    output << "assignment_epsilon="
           << options.place.assignment_recovery.assignment_epsilon << '\n';
    output << "assignment_batches="
           << result.assignment_recovery.assignment_batches << '\n';
    output << "assigned_nodes=" << result.assignment_recovery.assigned_nodes << '\n';
    output << "assignment_bids=" << result.assignment_recovery.assignment_bids << '\n';
    output << "assignment_final_hpwl="
           << result.assignment_recovery.final_hpwl << '\n';
    output << "assignment_final_overflow="
           << result.assignment_recovery.final_overflow << '\n';
    output << "assignment_seconds="
           << result.assignment_recovery.wall_seconds << '\n';
    output << "selected_iteration=" << result.selected_iteration << '\n';
    output << "feasible=" << (result.feasible ? "true" : "false") << '\n';
    output << "gp_hpwl=" << result.final_metrics.hpwl << '\n';
    output << "gp_overflow=" << result.final_metrics.overflow << '\n';
    output << "max_density=" << result.final_metrics.max_density << '\n';
    output << "density_energy=" << result.final_metrics.density_energy << '\n';
    output << "objective_evaluations=" << result.objective_evaluations << '\n';
    output << "gp_seconds=" << result.wall_seconds << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        if (!options.initial_placement.empty()) {
            std::string initial_name = options.initial_placement.string();
            std::transform(initial_name.begin(), initial_name.end(), initial_name.begin(),
                           [](unsigned char value) {
                               return static_cast<char>(std::tolower(value));
                           });
            if (initial_name.find("eplace") != std::string::npos) {
                throw std::invalid_argument(
                    "ePlace initial placements are disabled; use a non-ePlace "
                    "placement or omit --initial-placement for center initialization");
            }
        }
        ea::Database db = ea::read_bookshelf(options.benchmark);
        if (options.initial_placement.empty()) {
            ea::initialize_center_gaussian(db, options.seed, options.sigma_ratio);
        } else {
            ea::load_bookshelf_placement(db, options.initial_placement);
        }
        std::cout << "[Input] nodes=" << db.nodes.size()
                  << " movable=" << db.movable_ids.size()
                  << " fixed=" << db.fixed_ids.size()
                  << " nets=" << db.nets.size() << '\n';
        ea::PlaceResult result;
        if (!options.multilevel) {
            result = ea::global_place(db, options.place);
        } else {
            if (options.level_iterations <= 0) {
                throw std::invalid_argument("--level-iterations must be positive");
            }
            // This is a single staged pipeline: only the placement produced
            // by the preceding exact level is used by the next level. Every
            // level rebuilds fixed macro capacity and exact overlap from the
            // current rectangles; no occupancy interpolation is used.
            const std::vector<int> levels{64, 128, 256, 512};
            std::vector<ea::Real> staged_regional_prices;
            for (std::size_t level = 0; level < levels.size(); ++level) {
                ea::PlaceConfig stage = options.place;
                stage.bins_x = stage.bins_y = levels[level];
                stage.iterations = options.level_iterations;
                stage.output_dir = options.output /
                    ("level_" + std::to_string(levels[level]));
                stage.lambda.initial_effective = 0.0;
                // Density gradients change units with the bin pitch. Each
                // exact level therefore rebalances lambda from its own
                // gradients instead of inheriting an incompatible control.
                stage.lambda.initial_control = 1.0;
                if (stage.regional_price.enabled &&
                    !staged_regional_prices.empty()) {
                    stage.regional_price.initial_prices =
                        staged_regional_prices;
                }
                std::filesystem::create_directories(stage.output_dir);
                std::cout << "[Level] bins=" << levels[level]
                          << " iterations=" << stage.iterations << '\n';
                result = ea::global_place(db, stage);
                if (stage.regional_price.enabled && level + 1 < levels.size()) {
                    staged_regional_prices = ea::prolongate_bin_field(
                        result.final_regional_prices, levels[level],
                        levels[level], levels[level + 1], levels[level + 1]);
                }
                ea::write_bookshelf_placement(
                    db, stage.output_dir / "global.pl");
            }
        }
        ea::write_bookshelf_placement(db, options.output / "global.pl");
        write_summary(options, db, result);
        std::cout << "[Result] hpwl=" << result.final_metrics.hpwl
                  << " overflow=" << result.final_metrics.overflow
                  << " selected_iteration=" << result.selected_iteration
                  << " seconds=" << result.wall_seconds << '\n';
        return result.feasible ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
