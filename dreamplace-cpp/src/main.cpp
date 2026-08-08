#include "bookshelf.h"
#include "electric.h"
#include "legalizer.h"
#include "optimizer.h"
#include "wirelength.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;
using namespace dpcpp;

namespace {

struct Options {
    std::string benchmark;
    std::string output_dir;
    std::string initial_pl;
    GlobalPlaceConfig gp;
    LegalizeConfig legal;
    Real sigma_ratio = 0.001;
};

void usage() {
    std::cout
        << "Usage: dreamplace_cpp.exe <raw-benchmark-base> [options]\n"
        << "  Example: dreamplace_cpp.exe ../alg-electronic/ispd2005/adaptec1/adaptec1\n"
        << "Options:\n"
        << "  --iterations N       global-placement iterations (default 1000)\n"
        << "  --bins N             NxN electrostatic grid\n"
        << "  --bins-x N           horizontal bins\n"
        << "  --bins-y N           vertical bins\n"
        << "  --target-density X   density target (default 1.0)\n"
        << "  --stop-overflow X    feasible checkpoint threshold (default 0.07)\n"
        << "  --seed N             random seed (default 1000)\n"
        << "  --sigma-ratio X      center Gaussian sigma / region size (default 0.001)\n"
        << "  --initial-pl FILE    initialize movable cells from an explicit PL file\n"
        << "  --output DIR         output directory\n"
        << "  --no-fillers         disable filler cells\n"
        << "  --no-abacus          stop legalization after Greedy\n"
        << "  --no-detailed        disable legal detailed refinement\n"
        << "  --dreamplace-detailed  enable K-Reorder, global swap and independent-set matching\n"
        << "  --k-reorder-size N  K-Reorder window size (2..6, default 4)\n"
        << "  --dp-passes N       detailed-placement pass count\n"
        << "  --log-every N        metric print interval\n"
        << "  --exact-hpwl         exact nonsmooth HPWL objective (default optimizer Adam)\n"
        << "  --smooth-hpwl        weighted-average HPWL objective (default)\n"
        << "  --optimizer NAME     dreamplace|heavy-ball|adam|amsgrad|adagrad\n"
        << "  --lambda-policy NAME dreamplace|trajectory|ratio\n"
        << "  --hpwl-only N        exact-HPWL iterations before density activation\n"
        << "  --step-fraction X    nonsmooth density-stage step / region perimeter\n"
        << "  --hpwl-step-fraction X  nonsmooth HPWL-only step / region perimeter\n"
        << "  --lambda-interval N  trajectory/ratio update interval\n"
        << "  --trajectory-horizon N  iterations used by the overflow trajectory\n"
        << "  --hpwl-baseline X    HPWL reference for lambda feedback\n"
        << "  --overflow-baseline X  overflow reference for lambda feedback\n"
        << "  --lambda-min X       minimum dimensionless lambda control\n"
        << "  --lambda-max X       maximum dimensionless lambda control\n"
        << "  --bundle-hpwl        aggregate exact-HPWL cutting planes\n"
        << "  --bundle-size N      maximum retained HPWL cuts\n"
        << "  --bundle-groups N    experimental independent net groups (default 1)\n"
        << "  --bundle-group-size N  cuts retained per experimental net group\n"
        << "  --bundle-interval N  regular cut insertion interval\n"
        << "  --bundle-start-overflow X  enable bundle below this overflow\n"
        << "  --bundle-prox X      bundle proximal scaling\n"
        << "  --bundle-current-mix X  current-subgradient fraction in [0,1]\n"
        << "  --feasible-refinement  reset moments and refine along the 6.5%-7% boundary\n"
        << "  --refine-lower-overflow X  lower edge of feasible band (default 0.065)\n"
        << "  --refine-lr-scale X  post-feasible learning-rate multiplier (default 0.05)\n"
        << "  --refine-lambda-gain X  feasible-band lambda feedback gain (default 0.20)\n"
        << "  --snapshot-every N  save global-layout snapshots every N steps\n"
        << "  --snapshot-dir DIR  snapshot output directory (default output/snapshots)\n"
        << "Experimental opt-in controls (defaults preserve the original flow):\n"
        << "  --legal-checkpoints  rank feasible GP checkpoints by quick legal HPWL\n"
        << "  --legal-checkpoint-interval N  quick-legal evaluation interval\n"
        << "  --late-legal-projection  blend late GP states toward a legal projection\n"
        << "  --legal-projection-interval N  late projection interval\n"
        << "  --legal-projection-mix X  legal projection blend in [0,1]\n"
        << "  --legal-row-force X  late exact-HPWL row-proximity force\n"
        << "  --mixed-spectral-field  use mixed sine/cosine electric field\n"
        << "  --multilevel-density  progress from coarse to requested density grid\n"
        << "  --multilevel-min-bins N  initial grid for multilevel density\n"
        << "  --gradient-samples N  nearby exact-subgradient samples per update\n"
        << "  --gradient-sampling-interval N  gradient-sampling interval\n"
        << "  --gradient-sampling-radius X  sampling radius in placement units\n"
        << "  --serious-bundle  enable serious/null trial acceptance below bundle threshold\n"
        << "  --serious-step-ratio X  actual/predicted decrease threshold\n"
        << "  --bundle-overflow-tolerance X  trial-only overflow filter tolerance\n"
        << "  --progressive-legalization  enable the opt-in exact five-stage flow\n"
        << "  --progressive-density-iterations N  phase-2 iterations before macro ramp\n"
        << "  --progressive-overflow-lower X  lower density-band edge (default 0.07)\n"
        << "  --progressive-overflow-upper X  upper density-band edge (default 0.08)\n"
        << "  --progressive-hpwl-target X  HPWL budget used by band-dual control\n"
        << "  --progressive-dual-kp X  proportional log-lambda gain\n"
        << "  --progressive-dual-ki X  integral log-lambda gain\n"
        << "  --progressive-hpwl-gain X  HPWL budget release gain\n"
        << "  --progressive-dual-switch-overflow X  switch to band-dual control\n"
        << "  --progressive-obstacle-scale X  auto-scaled macro-force fraction\n"
        << "  --progressive-row-force X  final continuous row-force strength\n"
        << "  --progressive-segment-force X  final obstacle-free segment force\n"
        << "  --progressive-congestion-gain X  soft overfull-segment assignment cost\n"
        << "  --progressive-obstacle-iterations N  phase-3 force ramp length\n"
        << "  --progressive-filter-backtracks N  trial-filter backtracking limit\n"
        << "  --legal-refine-rounds N  outer legal detailed-placement rounds\n"
        << "  --cell-insertion-passes N  long-range legal insertion passes\n"
        << "  --cell-insertion-window N  maximum insertion rank distance\n"
        << "  --legal-projected-passes N  legal projected exact-HPWL passes\n"
        << "  --legal-projected-step-sites X  initial projected step in sites\n"
        << "  --legal-bundle-passes N  fixed-topology constrained bundle passes\n"
        << "  --legal-bundle-size N  retained cuts in the legal bundle\n"
        << "  --legal-bundle-step-sites X  initial legal-bundle step in sites\n"
        << "  --row-relegalization-passes N  adjacent-row relegalization passes\n"
        << "  --independent-set-size N  independent matching group size\n"
        << "  --hungarian-matching  use polynomial independent-set assignment\n";
}

std::string require_value(int& i, int argc, char** argv) {
    if (++i >= argc) throw std::runtime_error(std::string("missing value after ") + argv[i - 1]);
    return argv[i];
}

Options parse_options(int argc, char** argv) {
    if (argc < 2 || std::string(argv[1]) == "--help") {
        usage();
        std::exit(argc < 2 ? 1 : 0);
    }
    Options options;
    options.benchmark = argv[1];
    if (fs::path(options.benchmark).extension() == ".aux") {
        options.benchmark = fs::path(options.benchmark).replace_extension().string();
    }
    const std::string name = fs::path(options.benchmark).filename().string();
    if (name == "adaptec2") options.gp.bins_x = options.gp.bins_y = 1024;
    bool optimizer_explicit = false;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--iterations") options.gp.iterations = std::stoi(require_value(i, argc, argv));
        else if (arg == "--bins") {
            const int bins = std::stoi(require_value(i, argc, argv));
            options.gp.bins_x = options.gp.bins_y = bins;
        } else if (arg == "--bins-x") options.gp.bins_x = std::stoi(require_value(i, argc, argv));
        else if (arg == "--bins-y") options.gp.bins_y = std::stoi(require_value(i, argc, argv));
        else if (arg == "--target-density") options.gp.target_density = std::stod(require_value(i, argc, argv));
        else if (arg == "--stop-overflow") options.gp.stop_overflow = std::stod(require_value(i, argc, argv));
        else if (arg == "--seed") options.gp.seed = std::stoull(require_value(i, argc, argv));
        else if (arg == "--sigma-ratio") options.sigma_ratio = std::stod(require_value(i, argc, argv));
        else if (arg == "--initial-pl") options.initial_pl = require_value(i, argc, argv);
        else if (arg == "--output") options.output_dir = require_value(i, argc, argv);
        else if (arg == "--log-every") options.gp.log_every = std::stoi(require_value(i, argc, argv));
        else if (arg == "--no-fillers") options.gp.enable_fillers = false;
        else if (arg == "--no-abacus") options.legal.run_abacus = false;
        else if (arg == "--no-detailed") options.legal.run_detailed = false;
        else if (arg == "--dreamplace-detailed") options.legal.run_dreamplace_detailed = true;
        else if (arg == "--k-reorder-size") options.legal.k_reorder_size = std::stoi(require_value(i, argc, argv));
        else if (arg == "--dp-passes") {
            const int passes = std::stoi(require_value(i, argc, argv));
            options.legal.detailed_passes = passes;
            options.legal.global_swap_passes = passes;
        }
        else if (arg == "--exact-hpwl") options.gp.wirelength_model = WirelengthModel::ExactHpwl;
        else if (arg == "--smooth-hpwl") options.gp.wirelength_model = WirelengthModel::WeightedAverage;
        else if (arg == "--optimizer") {
            const std::string value = require_value(i, argc, argv);
            optimizer_explicit = true;
            if (value == "dreamplace") options.gp.optimizer = GlobalOptimizer::DreamplaceNesterov;
            else if (value == "heavy-ball") options.gp.optimizer = GlobalOptimizer::HeavyBall;
            else if (value == "adam") options.gp.optimizer = GlobalOptimizer::Adam;
            else if (value == "amsgrad") options.gp.optimizer = GlobalOptimizer::AMSGrad;
            else if (value == "adagrad") options.gp.optimizer = GlobalOptimizer::AdaGrad;
            else throw std::runtime_error("unknown optimizer: " + value);
        } else if (arg == "--lambda-policy") {
            const std::string value = require_value(i, argc, argv);
            if (value == "dreamplace") options.gp.lambda_policy = LambdaPolicy::Dreamplace;
            else if (value == "trajectory") options.gp.lambda_policy = LambdaPolicy::Trajectory;
            else if (value == "ratio") options.gp.lambda_policy = LambdaPolicy::Ratio;
            else throw std::runtime_error("unknown lambda policy: " + value);
        } else if (arg == "--hpwl-only") options.gp.hpwl_only_iterations = std::stoi(require_value(i, argc, argv));
        else if (arg == "--step-fraction") options.gp.nonsmooth_step_fraction = std::stod(require_value(i, argc, argv));
        else if (arg == "--hpwl-step-fraction") options.gp.nonsmooth_hpwl_step_fraction = std::stod(require_value(i, argc, argv));
        else if (arg == "--lambda-interval") options.gp.lambda_update_interval = std::stoi(require_value(i, argc, argv));
        else if (arg == "--trajectory-horizon") options.gp.lambda_trajectory_horizon = std::stoi(require_value(i, argc, argv));
        else if (arg == "--hpwl-baseline") options.gp.hpwl_baseline = std::stod(require_value(i, argc, argv));
        else if (arg == "--overflow-baseline") options.gp.overflow_baseline = std::stod(require_value(i, argc, argv));
        else if (arg == "--lambda-min") options.gp.lambda_control_min = std::stod(require_value(i, argc, argv));
        else if (arg == "--lambda-max") options.gp.lambda_control_max = std::stod(require_value(i, argc, argv));
        else if (arg == "--bundle-hpwl") options.gp.enable_bundle = true;
        else if (arg == "--bundle-size") options.gp.bundle_size = std::stoi(require_value(i, argc, argv));
        else if (arg == "--bundle-groups") options.gp.bundle_groups = std::stoi(require_value(i, argc, argv));
        else if (arg == "--bundle-group-size") options.gp.bundle_group_size = std::stoi(require_value(i, argc, argv));
        else if (arg == "--bundle-interval") options.gp.bundle_interval = std::stoi(require_value(i, argc, argv));
        else if (arg == "--bundle-start-overflow") options.gp.bundle_start_overflow = std::stod(require_value(i, argc, argv));
        else if (arg == "--bundle-prox") options.gp.bundle_prox_scale = std::stod(require_value(i, argc, argv));
        else if (arg == "--bundle-current-mix") options.gp.bundle_current_mix = std::stod(require_value(i, argc, argv));
        else if (arg == "--feasible-refinement") options.gp.feasible_refinement = true;
        else if (arg == "--refine-lower-overflow") options.gp.refinement_lower_overflow = std::stod(require_value(i, argc, argv));
        else if (arg == "--refine-lr-scale") options.gp.refinement_learning_rate_scale = std::stod(require_value(i, argc, argv));
        else if (arg == "--refine-lambda-gain") options.gp.refinement_lambda_gain = std::stod(require_value(i, argc, argv));
        else if (arg == "--snapshot-every") options.gp.snapshot_every = std::stoi(require_value(i, argc, argv));
        else if (arg == "--snapshot-dir") options.gp.snapshot_dir = require_value(i, argc, argv);
        else if (arg == "--legal-checkpoints") options.gp.legal_checkpoint_selection = true;
        else if (arg == "--legal-checkpoint-interval") options.gp.legal_checkpoint_interval = std::stoi(require_value(i, argc, argv));
        else if (arg == "--late-legal-projection") options.gp.late_legal_projection = true;
        else if (arg == "--legal-projection-interval") options.gp.legal_projection_interval = std::stoi(require_value(i, argc, argv));
        else if (arg == "--legal-projection-mix") options.gp.legal_projection_mix = std::stod(require_value(i, argc, argv));
        else if (arg == "--legal-row-force") options.gp.legal_row_force = std::stod(require_value(i, argc, argv));
        else if (arg == "--mixed-spectral-field") options.gp.mixed_spectral_field = true;
        else if (arg == "--multilevel-density") options.gp.multilevel_density = true;
        else if (arg == "--multilevel-min-bins") options.gp.multilevel_min_bins = std::stoi(require_value(i, argc, argv));
        else if (arg == "--gradient-samples") options.gp.gradient_sampling_samples = std::stoi(require_value(i, argc, argv));
        else if (arg == "--gradient-sampling-interval") options.gp.gradient_sampling_interval = std::stoi(require_value(i, argc, argv));
        else if (arg == "--gradient-sampling-radius") options.gp.gradient_sampling_radius = std::stod(require_value(i, argc, argv));
        else if (arg == "--serious-bundle") options.gp.serious_bundle = true;
        else if (arg == "--serious-step-ratio") options.gp.serious_step_ratio = std::stod(require_value(i, argc, argv));
        else if (arg == "--bundle-overflow-tolerance") options.gp.bundle_overflow_tolerance = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-legalization") options.gp.progressive_legalization = true;
        else if (arg == "--progressive-density-iterations") options.gp.progressive_density_iterations = std::stoi(require_value(i, argc, argv));
        else if (arg == "--progressive-overflow-lower") options.gp.progressive_overflow_lower = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-overflow-upper") options.gp.progressive_overflow_upper = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-hpwl-target") options.gp.progressive_hpwl_target = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-dual-kp") options.gp.progressive_dual_kp = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-dual-ki") options.gp.progressive_dual_ki = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-hpwl-gain") options.gp.progressive_hpwl_gain = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-dual-switch-overflow") options.gp.progressive_dual_switch_overflow = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-obstacle-scale") options.gp.progressive_obstacle_scale = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-row-force") options.gp.progressive_row_force = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-segment-force") options.gp.progressive_segment_force = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-congestion-gain") options.gp.progressive_congestion_gain = std::stod(require_value(i, argc, argv));
        else if (arg == "--progressive-obstacle-iterations") options.gp.progressive_obstacle_iterations = std::stoi(require_value(i, argc, argv));
        else if (arg == "--progressive-filter-backtracks") options.gp.progressive_filter_backtracks = std::stoi(require_value(i, argc, argv));
        else if (arg == "--legal-refine-rounds") options.legal.detailed_outer_rounds = std::stoi(require_value(i, argc, argv));
        else if (arg == "--cell-insertion-passes") options.legal.cell_insertion_passes = std::stoi(require_value(i, argc, argv));
        else if (arg == "--cell-insertion-window") options.legal.cell_insertion_window = std::stoi(require_value(i, argc, argv));
        else if (arg == "--legal-projected-passes") options.legal.projected_subgradient_passes = std::stoi(require_value(i, argc, argv));
        else if (arg == "--legal-projected-step-sites") options.legal.projected_step_sites = std::stod(require_value(i, argc, argv));
        else if (arg == "--legal-bundle-passes") options.legal.constrained_bundle_passes = std::stoi(require_value(i, argc, argv));
        else if (arg == "--legal-bundle-size") options.legal.constrained_bundle_size = std::stoi(require_value(i, argc, argv));
        else if (arg == "--legal-bundle-step-sites") options.legal.constrained_bundle_step_sites = std::stod(require_value(i, argc, argv));
        else if (arg == "--row-relegalization-passes") options.legal.row_relegalization_passes = std::stoi(require_value(i, argc, argv));
        else if (arg == "--independent-set-size") options.legal.independent_set_size = std::stoi(require_value(i, argc, argv));
        else if (arg == "--hungarian-matching") options.legal.use_hungarian_matching = true;
        else throw std::runtime_error("unknown option: " + arg);
    }
    if (!optimizer_explicit) {
        options.gp.optimizer = options.gp.progressive_legalization
            ? GlobalOptimizer::AMSGrad
            : options.gp.wirelength_model == WirelengthModel::ExactHpwl
            ? GlobalOptimizer::Adam : GlobalOptimizer::DreamplaceNesterov;
    }
    if (options.gp.progressive_legalization) {
        options.gp.wirelength_model = WirelengthModel::ExactHpwl;
        options.gp.lambda_policy = LambdaPolicy::BandDual;
    }
    if (options.gp.iterations <= 0 || options.gp.bins_x <= 0 || options.gp.bins_y <= 0)
        throw std::runtime_error("iterations and bin counts must be positive");
    if (options.gp.hpwl_only_iterations < 0 || options.gp.lambda_update_interval < 0 ||
        options.gp.bundle_size <= 0 || options.gp.bundle_interval <= 0 ||
        options.gp.bundle_groups <= 0 || options.gp.bundle_group_size <= 0)
        throw std::runtime_error("iteration intervals and bundle size are invalid");
    if (options.gp.nonsmooth_step_fraction <= 0.0 ||
        options.gp.nonsmooth_hpwl_step_fraction <= 0.0)
        throw std::runtime_error("nonsmooth step fractions must be positive");
    if (options.gp.bundle_current_mix < 0.0 || options.gp.bundle_current_mix > 1.0)
        throw std::runtime_error("bundle current mix must lie in [0,1]");
    if (options.gp.lambda_control_min <= 0.0 ||
        options.gp.lambda_control_max < options.gp.lambda_control_min)
        throw std::runtime_error("lambda control bounds are invalid");
    if (options.gp.refinement_lower_overflow < 0.0 ||
        options.gp.refinement_lower_overflow >= options.gp.stop_overflow ||
        options.gp.refinement_learning_rate_scale <= 0.0 ||
        options.gp.refinement_lambda_gain <= 0.0)
        throw std::runtime_error("feasible-refinement parameters are invalid");
    if (options.gp.snapshot_every < 0)
        throw std::runtime_error("snapshot interval must be non-negative");
    if (options.gp.legal_checkpoint_interval <= 0 ||
        options.gp.legal_projection_interval <= 0 ||
        options.gp.gradient_sampling_interval <= 0 ||
        options.gp.gradient_sampling_samples < 0 ||
        options.gp.gradient_sampling_radius < 0.0 ||
        options.gp.multilevel_min_bins <= 0)
        throw std::runtime_error("experimental GP intervals and sizes are invalid");
    if (options.gp.legal_projection_mix < 0.0 || options.gp.legal_projection_mix > 1.0 ||
        options.gp.legal_projection_max_hpwl_ratio < 1.0 ||
        options.gp.legal_row_force < 0.0 || options.gp.serious_step_ratio < 0.0 ||
        options.gp.serious_step_ratio > 1.0 || options.gp.bundle_overflow_tolerance < 0.0)
        throw std::runtime_error("experimental GP gains are invalid");
    if (options.gp.progressive_density_iterations < 0 ||
        options.gp.progressive_overflow_lower < 0.0 ||
        options.gp.progressive_overflow_upper <= options.gp.progressive_overflow_lower ||
        options.gp.progressive_hpwl_target <= 0.0 ||
        options.gp.progressive_dual_kp < 0.0 ||
        options.gp.progressive_dual_ki < 0.0 ||
        options.gp.progressive_hpwl_gain < 0.0 ||
        options.gp.progressive_dual_switch_overflow <= options.gp.progressive_overflow_upper ||
        options.gp.progressive_obstacle_scale < 0.0 ||
        options.gp.progressive_row_force < 0.0 ||
        options.gp.progressive_segment_force < 0.0 ||
        options.gp.progressive_congestion_gain < 0.0 ||
        options.gp.progressive_obstacle_iterations <= 0 ||
        options.gp.progressive_filter_backtracks < 0)
        throw std::runtime_error("progressive-legalization parameters are invalid");
    if (options.legal.k_reorder_size < 2 || options.legal.k_reorder_size > 6 ||
        options.legal.detailed_passes < 0 || options.legal.global_swap_passes < 0 ||
        options.legal.detailed_outer_rounds <= 0 ||
        options.legal.cell_insertion_passes < 0 || options.legal.cell_insertion_window <= 0 ||
        options.legal.projected_subgradient_passes < 0 ||
        options.legal.projected_step_sites <= 0.0 ||
        options.legal.constrained_bundle_passes < 0 ||
        options.legal.constrained_bundle_size <= 0 ||
        options.legal.constrained_bundle_step_sites <= 0.0 ||
        options.legal.row_relegalization_passes < 0 ||
        options.legal.independent_set_size < 2 || options.legal.independent_set_size > 32)
        throw std::runtime_error("detailed-placement parameters are invalid");
    if (options.output_dir.empty()) options.output_dir = (fs::path("output") / name).string();
    if (options.gp.snapshot_dir.empty()) options.gp.snapshot_dir =
        (fs::path(options.output_dir) / "snapshots").string();
    return options;
}

void write_summary(const fs::path& path, const Database& db,
                   const GlobalPlaceResult& gp, const LegalizeResult& legal,
                   const Options& options) {
    std::ofstream out(path);
    out << std::setprecision(12);
    out << "benchmark=" << fs::path(options.benchmark).filename().string() << '\n';
    out << "raw_pl=" << db.raw_pl_path << '\n';
    out << "initialization=" << (options.initial_pl.empty() ? "center_gaussian" : "external_pl") << '\n';
    if (!options.initial_pl.empty()) {
        out << "initial_pl=" << fs::absolute(options.initial_pl).lexically_normal().string() << '\n';
    }
    out << "seed=" << options.gp.seed << '\n';
    out << "sigma_ratio=" << options.sigma_ratio << '\n';
    out << "bins=" << options.gp.bins_x << 'x' << options.gp.bins_y << '\n';
    out << "iterations=" << options.gp.iterations << '\n';
    out << "wirelength_model=" << wirelength_model_name(options.gp.wirelength_model) << '\n';
    out << "optimizer=" << global_optimizer_name(options.gp.optimizer) << '\n';
    out << "lambda_policy=" << lambda_policy_name(options.gp.lambda_policy) << '\n';
    out << "lambda_control_min=" << options.gp.lambda_control_min << '\n';
    out << "lambda_control_max=" << options.gp.lambda_control_max << '\n';
    out << "hpwl_only_iterations=" << options.gp.hpwl_only_iterations << '\n';
    out << "bundle_hpwl=" << (options.gp.enable_bundle ? "true" : "false") << '\n';
    out << "bundle_groups=" << options.gp.bundle_groups << '\n';
    out << "bundle_group_size=" << options.gp.bundle_group_size << '\n';
    out << "bundle_overflow_tolerance=" << options.gp.bundle_overflow_tolerance << '\n';
    out << "feasible_refinement=" << (options.gp.feasible_refinement ? "true" : "false") << '\n';
    out << "legal_checkpoint_selection=" << (options.gp.legal_checkpoint_selection ? "true" : "false") << '\n';
    out << "late_legal_projection=" << (options.gp.late_legal_projection ? "true" : "false") << '\n';
    out << "mixed_spectral_field=" << (options.gp.mixed_spectral_field ? "true" : "false") << '\n';
    out << "multilevel_density=" << (options.gp.multilevel_density ? "true" : "false") << '\n';
    out << "gradient_sampling_samples=" << options.gp.gradient_sampling_samples << '\n';
    out << "serious_bundle=" << (options.gp.serious_bundle ? "true" : "false") << '\n';
    out << "progressive_legalization=" << (options.gp.progressive_legalization ? "true" : "false") << '\n';
    out << "progressive_density_iterations=" << options.gp.progressive_density_iterations << '\n';
    out << "progressive_overflow_band=" << options.gp.progressive_overflow_lower
        << ',' << options.gp.progressive_overflow_upper << '\n';
    out << "progressive_hpwl_target=" << options.gp.progressive_hpwl_target << '\n';
    out << "progressive_dual_switch_overflow=" << options.gp.progressive_dual_switch_overflow << '\n';
    out << "progressive_obstacle_scale=" << options.gp.progressive_obstacle_scale << '\n';
    out << "progressive_row_force=" << options.gp.progressive_row_force << '\n';
    out << "progressive_segment_force=" << options.gp.progressive_segment_force << '\n';
    out << "progressive_congestion_gain=" << options.gp.progressive_congestion_gain << '\n';
    out << "progressive_obstacle_iterations=" << options.gp.progressive_obstacle_iterations << '\n';
    out << "legal_refine_rounds=" << options.legal.detailed_outer_rounds << '\n';
    out << "cell_insertion_passes=" << options.legal.cell_insertion_passes << '\n';
    out << "projected_subgradient_passes=" << options.legal.projected_subgradient_passes << '\n';
    out << "constrained_bundle_passes=" << options.legal.constrained_bundle_passes << '\n';
    out << "constrained_bundle_size=" << options.legal.constrained_bundle_size << '\n';
    out << "row_relegalization_passes=" << options.legal.row_relegalization_passes << '\n';
    out << "independent_set_size=" << options.legal.independent_set_size << '\n';
    out << "gp_hpwl=" << gp.final_metrics.exact_hpwl << '\n';
    out << "gp_overflow=" << gp.final_metrics.overflow << '\n';
    out << "gp_seconds=" << gp.wall_time_seconds << '\n';
    out << "selected_legal_hpwl=" << gp.selected_legal_hpwl << '\n';
    out << "selected_legal_iteration=" << gp.selected_legal_iteration << '\n';
    out << "prelegal_hpwl=" << legal.hpwl_before << '\n';
    out << "greedy_hpwl=" << legal.hpwl_after_greedy << '\n';
    out << "abacus_hpwl=" << legal.hpwl_after_abacus << '\n';
    out << "k_reorder_hpwl=" << legal.hpwl_after_k_reorder << '\n';
    out << "global_swap_hpwl=" << legal.hpwl_after_global_swap << '\n';
    out << "independent_set_hpwl=" << legal.hpwl_after_independent_set << '\n';
    out << "constrained_bundle_hpwl=" << legal.hpwl_after_constrained_bundle << '\n';
    out << "detailed_hpwl=" << legal.hpwl_after_detailed << '\n';
    out << "internal_legal=" << (legal.legality.legal ? "true" : "false") << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        fs::create_directories(options.output_dir);
        Database db = read_bookshelf(options.benchmark);
        if (options.initial_pl.empty()) {
            center_gaussian_initialize(db, options.gp.seed, options.sigma_ratio);
        } else {
            load_movable_placement(db, options.initial_pl);
        }

        std::vector<Filler> fillers;
        if (options.gp.enable_fillers) {
            fillers = initialize_fillers(db, options.gp.target_density, options.gp.seed);
        }
        GlobalPlaceConfig gp_config = options.gp;
        if (!options.initial_pl.empty()) gp_config.gp_noise_ratio = 0.0;
        gp_config.metrics_path = (fs::path(options.output_dir) / "global_metrics.csv").string();
        if (gp_config.snapshot_every > 0) fs::create_directories(gp_config.snapshot_dir);
        const GlobalPlaceResult gp = global_place(db, fillers, gp_config);
        write_bookshelf_pl(db, (fs::path(options.output_dir) / "global.pl").string());
        if (gp_config.snapshot_every > 0) {
            write_bookshelf_pl(db, (fs::path(gp_config.snapshot_dir) / "global_selected.pl").string());
        }

        LegalizeConfig legal_config = options.legal;
        legal_config.snapshot_dir = gp_config.snapshot_every > 0 ? gp_config.snapshot_dir : "";
        const LegalizeResult legal = legalize_and_refine(db, legal_config);
        const fs::path final_path = fs::path(options.output_dir) / "final.pl";
        write_bookshelf_pl(db, final_path.string());
        write_summary(fs::path(options.output_dir) / "summary.txt", db, gp, legal, options);

        std::cout << "[Result] final_pl=" << fs::absolute(final_path).string() << '\n'
                  << "[Result] GP HPWL=" << gp.final_metrics.exact_hpwl
                  << " overflow=" << gp.final_metrics.overflow << '\n'
                  << "[Result] final legal HPWL=" << legal.hpwl_after_detailed
                  << " internal_legal=" << (legal.legality.legal ? "yes" : "no") << '\n';
        return legal.legality.legal ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
