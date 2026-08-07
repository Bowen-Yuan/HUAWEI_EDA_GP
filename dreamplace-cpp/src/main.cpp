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
        << "  --bundle-interval N  regular cut insertion interval\n"
        << "  --bundle-start-overflow X  enable bundle below this overflow\n"
        << "  --bundle-prox X      bundle proximal scaling\n"
        << "  --bundle-current-mix X  current-subgradient fraction in [0,1]\n"
        << "  --feasible-refinement  reset moments and refine along the 6.5%-7% boundary\n"
        << "  --refine-lower-overflow X  lower edge of feasible band (default 0.065)\n"
        << "  --refine-lr-scale X  post-feasible learning-rate multiplier (default 0.05)\n"
        << "  --refine-lambda-gain X  feasible-band lambda feedback gain (default 0.20)\n"
        << "  --snapshot-every N  save exact global-layout snapshots every N steps\n"
        << "  --snapshot-dir DIR  snapshot output directory (default output/snapshots)\n";
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
        else throw std::runtime_error("unknown option: " + arg);
    }
    if (!optimizer_explicit) {
        options.gp.optimizer = options.gp.wirelength_model == WirelengthModel::ExactHpwl
            ? GlobalOptimizer::Adam : GlobalOptimizer::DreamplaceNesterov;
    }
    if (options.gp.iterations <= 0 || options.gp.bins_x <= 0 || options.gp.bins_y <= 0)
        throw std::runtime_error("iterations and bin counts must be positive");
    if (options.gp.hpwl_only_iterations < 0 || options.gp.lambda_update_interval < 0 ||
        options.gp.bundle_size <= 0 || options.gp.bundle_interval <= 0)
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
    if (options.legal.k_reorder_size < 2 || options.legal.k_reorder_size > 6 ||
        options.legal.detailed_passes < 0 || options.legal.global_swap_passes < 0)
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
    out << "initialization=center_gaussian\n";
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
    out << "feasible_refinement=" << (options.gp.feasible_refinement ? "true" : "false") << '\n';
    out << "gp_hpwl=" << gp.final_metrics.exact_hpwl << '\n';
    out << "gp_overflow=" << gp.final_metrics.overflow << '\n';
    out << "gp_seconds=" << gp.wall_time_seconds << '\n';
    out << "prelegal_hpwl=" << legal.hpwl_before << '\n';
    out << "greedy_hpwl=" << legal.hpwl_after_greedy << '\n';
    out << "abacus_hpwl=" << legal.hpwl_after_abacus << '\n';
    out << "k_reorder_hpwl=" << legal.hpwl_after_k_reorder << '\n';
    out << "global_swap_hpwl=" << legal.hpwl_after_global_swap << '\n';
    out << "independent_set_hpwl=" << legal.hpwl_after_independent_set << '\n';
    out << "detailed_hpwl=" << legal.hpwl_after_detailed << '\n';
    out << "internal_legal=" << (legal.legality.legal ? "true" : "false") << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        fs::create_directories(options.output_dir);
        Database db = read_bookshelf(options.benchmark);
        center_gaussian_initialize(db, options.gp.seed, options.sigma_ratio);

        std::vector<Filler> fillers;
        if (options.gp.enable_fillers) {
            fillers = initialize_fillers(db, options.gp.target_density, options.gp.seed);
        }
        GlobalPlaceConfig gp_config = options.gp;
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
