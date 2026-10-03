#include "types.h"
#include "initial_placement.h"
#include "visualization.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <ctime>

// Forward declarations
bool parse_ispd2005(const std::string& bench_name,
                     std::vector<Cell>& cells,
                     std::vector<Net>& nets,
                     float_t& chip_xl, float_t& chip_yl,
                     float_t& chip_xh, float_t& chip_yh,
                     int& num_terminals);

PlacementResult solve_placement(std::vector<Cell>& cells,
                                 std::vector<Net>& nets,
                                 float_t chip_xl, float_t chip_yl,
                                 float_t chip_xh, float_t chip_yh,
                                 const SolverConfig& cfg);

struct CWTRConfig;
PlacementResult solve_cwtr(std::vector<Cell>& cells,
                            std::vector<Net>& nets,
                            float_t chip_xl, float_t chip_yl,
                            float_t chip_xh, float_t chip_yh,
                            const CWTRConfig& cfg);

PlacementResult solve_spread_first(std::vector<Cell>& cells,
                                    std::vector<Net>& nets,
                                    float_t chip_xl, float_t chip_yl,
                                    float_t chip_xh, float_t chip_yh,
                                    float_t target_avg_overflow = 0.01f);

PlacementResult solve_diffusion_spread(std::vector<Cell>& cells,
                                        std::vector<Net>& nets,
                                        float_t chip_xl, float_t chip_yl,
                                        float_t chip_xh, float_t chip_yh,
                                        float_t target_avg_overflow = 0.01f);

PlacementResult solve_uniform_start(std::vector<Cell>& cells,
                                     std::vector<Net>& nets,
                                     float_t chip_xl, float_t chip_yl,
                                     float_t chip_xh, float_t chip_yh);

float_t compute_hpwl_only(const std::vector<Cell>&, const std::vector<Net>&);

void set_short_edge_only(bool v);
void legalize_placement(std::vector<Cell>& cells,
                         float_t chip_xl, float_t chip_yl,
                         float_t chip_xh, float_t chip_yh,
                         const std::string& scl_file);

void write_placement(const std::string& out_file,
                     const std::vector<Cell>& cells,
                     float_t chip_xl, float_t chip_yl,
                     float_t chip_xh, float_t chip_yh);

int main(int argc, char* argv[]) {
    std::setvbuf(stdout, nullptr, _IOLBF, 0);
    if (argc < 2) {
        printf("Usage: %s <benchmark_name> [--cwtr|--spread|--diffuse|--uniform] [--init <strategy>] [--bo]\n", argv[0]);
        printf("  default   : Phase A subgradient only\n");
        printf("  --cwtr    : Phase A + CWTR refinement\n");
        printf("  --uniform : Uniform start + HPWL optimization\n");
        printf("  --spread  : Aggressive spread + CWTR recovery\n");
        printf("  --diffuse : Diffusion potential spreading\n");
        printf("  --bo      : Enable Bayesian Optimization for λ selection\n");
        printf("\nInitial placement strategies (--init <name>):\n");
        printf("  eplace    : ePlace QP (default, from .eplace-ip.pl)\n");
        printf("  uniform   : Uniform sqrt(N) grid\n");
        printf("  random    : Uniform random within chip area\n");
        printf("  row       : Pack cells along .scl rows\n");
        printf("  scaled    : ePlace coords expanded 1.3x from center\n");
        printf("  gaussian  : Gaussian around chip center\n");
        printf("  quadrant  : Place cells near terminal centroid\n");
        printf("\nPhase iteration controls:\n");
        printf("  --p0-iters N  : Phase 0 HPWL-only iterations (default 150, 0=skip)\n");
        printf("  --p1-iters N  : Phase 1 coarse-spread iterations (default 100)\n");
        printf("  --p2-iters N  : Phase 2 medium-balance iterations (default 150)\n");
        printf("  --max-iters N : total NSP iteration budget (default 800)\n");
        printf("  --p1-lambda X : Phase 1 initial density penalty λ (default 0.005)\n");
        printf("  --lambda-interval N : baseline controller interval (default 50)\n");
        printf("  --lambda-up X : controller increase factor (default 1.5)\n");
        printf("  --lambda-down X : controller decrease factor (default 0.5)\n");
        printf("  --lambda-hysteresis X : ratio deadband (default 0.05)\n");
        printf("  --lambda-max X : maximum control lambda (default 20)\n");
        printf("  --overflow-baseline X : overflow controller baseline (default 0.06)\n");
        printf("  --overflow-limit X : best-state feasibility limit (default: baseline)\n");
        printf("  --target-density X : override automatic feasible density\n");
        printf("  --adam : use Adam updates in density phases\n");
        printf("  HPWL gradients are always exact max-minus-min subgradients.\n");
        printf("  --lambda-guard : overflow guard with HPWL recovery pulses\n");
        printf("  --trajectory-lambda : use target-trajectory PI-D lambda feedback\n");
        printf("  --baseline-lambda : use continuous baseline feedback (default)\n");
        printf("  --lambda-horizon N : fine steps to reach target overflow (default 1150)\n");
        printf("  --trajectory-target X : trajectory terminal overflow (default 0.05998)\n");
        printf("  --lambda-guard-low X : baseline recovery target (default 0.0598)\n");
        printf("  --lambda-guard-floor X : minimum spread pulse (default 2)\n");
        printf("  --lambda-guard-interval N : fine-grid guard interval (default 25)\n");
        printf("  --fine-momentum X : fine-grid heavy-ball momentum (default 0.70)\n");
        printf("  --fine-cooldown-floor X : final LR multiplier (default 0.15)\n");
        printf("  --projected-hpwl : remove HPWL components that increase electric energy near feasibility\n");
        printf("  --adaptive-restart : reset local heavy-ball momentum on gradient conflict\n");
        printf("  --projection-margin X : project only within 0.06+-X overflow (default 0.002)\n");
        printf("  --bundle-hpwl : use a proximal cutting-plane bundle in fine refinement\n");
        printf("  --bundle-size N : retained exact-HPWL cuts (default 4)\n");
        printf("  --bundle-interval N : fine steps between new cuts (default 5)\n");
        printf("  --bundle-adaptive-cuts : insert cuts on exact-subgradient direction changes\n");
        printf("  --bundle-min-interval N : minimum adaptive cut age (default 2)\n");
        printf("  --bundle-cosine X : adaptive insertion cosine threshold (default 0.90)\n");
        printf("  --bundle-prox X : bundle proximal scale (default 1.0)\n");
        printf("  --bundle-current-mix X : exact current-cut fraction (default 0.35)\n");
        printf("  --bundle-extrema-mix : adapt current-cut fraction from subgradient cosine\n");
        printf("  --bundle-mix-gain X : cosine feedback gain (default 0.30)\n");
        printf("  --bundle-mix-decay N : linearly remove cosine feedback over N fine steps\n");
        printf("  --short-edge  : Zero out less-effective density gradient direction\n");
        printf("  --visualize [DIR] : export optional snapshots (default visualizations)\n");
        printf("  --visualize-every N : snapshot interval (default 10)\n");
        printf("\nCurrent default strategy: 3000 iters, p1-lambda=0.02, "
               "heavy-ball, exact HPWL, continuous lambda target=0.0598.\n");
        return 1;
    }

    std::string bench_name = argv[1];
    bool use_cwtr = false, use_spread = false, use_diffuse = false;
    bool use_uniform = false;
    bool use_bo = false;
    std::string init_strategy = "eplace";  // default
    int p0_iters = 150;          // Phase 0: HPWL-only
    int p1_iters = 100;          // Phase 1: coarse spread
    int p2_iters = 150;          // Phase 2: medium balance
    int max_iters = 3000;        // Total NSP budget
    float_t p1_lambda = 0.02f;   // Phase 1 initial density weight
    int lambda_interval = 50;
    float_t lambda_up = 1.15f;
    float_t lambda_down = 0.5f;
    float_t lambda_hysteresis = 0.05f;
    float_t lambda_max = 100.0f;
    float_t overflow_baseline = 0.06f;
    float_t overflow_limit = -1.0f;
    float_t target_density_override = -1.0f;
    bool use_adam = false;
    bool use_lambda_guard = true;
    bool use_trajectory_lambda = false;
    float_t lambda_guard_low = 0.0598f;
    float_t lambda_guard_floor = 0.1f;
    int lambda_guard_interval = 25;
    int lambda_trajectory_horizon = 1150;
    float_t lambda_trajectory_target = 0.05998f;
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
    bool use_short_edge = false; // Short-edge-only density gradient
    VisualizationConfig visualization;

    // Parse flags
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--cwtr") use_cwtr = true;
        else if (arg == "--spread") use_spread = true;
        else if (arg == "--diffuse") use_diffuse = true;
        else if (arg == "--uniform") use_uniform = true;
        else if (arg == "--bo") use_bo = true;
        else if (arg == "--init" && i + 1 < argc) {
            init_strategy = argv[++i];
        }
        else if (arg == "--p0-iters" && i + 1 < argc) {
            p0_iters = std::atoi(argv[++i]);
        }
        else if (arg == "--p1-iters" && i + 1 < argc) {
            p1_iters = std::atoi(argv[++i]);
        }
        else if (arg == "--p2-iters" && i + 1 < argc) {
            p2_iters = std::atoi(argv[++i]);
        }
        else if (arg == "--max-iters" && i + 1 < argc) {
            max_iters = std::max(0, std::atoi(argv[++i]));
        }
        else if (arg == "--p1-lambda" && i + 1 < argc) {
            p1_lambda = std::atof(argv[++i]);
        }
        else if (arg == "--lambda-interval" && i + 1 < argc) {
            lambda_interval = std::max(1, std::atoi(argv[++i]));
        }
        else if (arg == "--lambda-up" && i + 1 < argc) {
            lambda_up = std::max(1.0f, (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--lambda-down" && i + 1 < argc) {
            lambda_down = std::max(0.0f, std::min(1.0f, (float_t)std::atof(argv[++i])));
        }
        else if (arg == "--lambda-hysteresis" && i + 1 < argc) {
            lambda_hysteresis = std::max(0.0f, (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--lambda-max" && i + 1 < argc) {
            lambda_max = std::max(1.0e-4f, (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--overflow-baseline" && i + 1 < argc) {
            overflow_baseline = std::max(1.0e-4f,
                (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--overflow-limit" && i + 1 < argc) {
            overflow_limit = std::max(1.0e-4f,
                (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--target-density" && i + 1 < argc) {
            target_density_override = std::max(1.0e-4f,
                (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--adam") {
            use_adam = true;
        }
        else if (arg == "--heavy-ball") {
            use_adam = false;
        }
        else if (arg == "--lambda-guard") {
            use_lambda_guard = true;
        }
        else if (arg == "--no-lambda-guard") {
            use_lambda_guard = false;
        }
        else if (arg == "--trajectory-lambda") {
            use_trajectory_lambda = true;
        }
        else if (arg == "--baseline-lambda") {
            use_trajectory_lambda = false;
        }
        else if (arg == "--lambda-horizon" && i + 1 < argc) {
            lambda_trajectory_horizon = std::max(1, std::atoi(argv[++i]));
        }
        else if (arg == "--trajectory-target" && i + 1 < argc) {
            lambda_trajectory_target = std::max(1.0e-4f,
                (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--lambda-guard-low" && i + 1 < argc) {
            lambda_guard_low = std::max(0.0f, (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--lambda-guard-floor" && i + 1 < argc) {
            lambda_guard_floor = std::max(1.0e-4f,
                (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--lambda-guard-interval" && i + 1 < argc) {
            lambda_guard_interval = std::max(1, std::atoi(argv[++i]));
        }
        else if (arg == "--fine-momentum" && i + 1 < argc) {
            fine_momentum = std::max(0.0f, std::min(0.99f,
                (float_t)std::atof(argv[++i])));
        }
        else if (arg == "--fine-cooldown-floor" && i + 1 < argc) {
            fine_cooldown_floor = std::max(0.01f, std::min(1.0f,
                (float_t)std::atof(argv[++i])));
        }
        else if (arg == "--projected-hpwl") {
            use_projected_hpwl = true;
        }
        else if (arg == "--adaptive-restart") {
            use_adaptive_restart = true;
        }
        else if (arg == "--projection-margin" && i + 1 < argc) {
            projection_overflow_margin = std::max(0.0f,
                (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--bundle-hpwl") {
            use_hpwl_bundle = true;
        }
        else if (arg == "--bundle-size" && i + 1 < argc) {
            hpwl_bundle_size = std::max(1, std::min(12, std::atoi(argv[++i])));
        }
        else if (arg == "--bundle-interval" && i + 1 < argc) {
            hpwl_bundle_interval = std::max(1, std::atoi(argv[++i]));
        }
        else if (arg == "--bundle-adaptive-cuts") {
            hpwl_bundle_adaptive_cuts = true;
        }
        else if (arg == "--bundle-min-interval" && i + 1 < argc) {
            hpwl_bundle_min_interval = std::max(1, std::atoi(argv[++i]));
        }
        else if (arg == "--bundle-cosine" && i + 1 < argc) {
            hpwl_bundle_cosine_trigger = std::max(-1.0f, std::min(1.0f,
                (float_t)std::atof(argv[++i])));
        }
        else if (arg == "--bundle-prox" && i + 1 < argc) {
            hpwl_bundle_prox_scale = std::max(0.01f,
                (float_t)std::atof(argv[++i]));
        }
        else if (arg == "--bundle-current-mix" && i + 1 < argc) {
            hpwl_bundle_current_mix = std::max(0.0f, std::min(1.0f,
                (float_t)std::atof(argv[++i])));
        }
        else if (arg == "--bundle-extrema-mix") {
            hpwl_bundle_extrema_mix = true;
        }
        else if (arg == "--bundle-mix-gain" && i + 1 < argc) {
            hpwl_bundle_mix_gain = std::max(0.0f, std::min(1.0f,
                (float_t)std::atof(argv[++i])));
        }
        else if (arg == "--bundle-mix-decay" && i + 1 < argc) {
            hpwl_bundle_mix_decay_steps = std::max(0, std::atoi(argv[++i]));
        }
        else if (arg == "--short-edge") {
            use_short_edge = true;
        }
        else if (arg == "--visualize") {
            visualization.enabled = true;
            visualization.output_dir = "visualizations";
            if (i + 1 < argc && argv[i + 1][0] != '-') visualization.output_dir = argv[++i];
        }
        else if (arg == "--visualize-every" && i + 1 < argc) {
            visualization.every = std::max(1, std::atoi(argv[++i]));
        }
    }
    if (visualization.enabled && visualization.every == 0) visualization.every = 10;

    printf("==========================================================\n");
    printf("  NSP: Non-Smooth Global Placer\n");
    printf("  Benchmark: %s", bench_name.c_str());
    if (use_uniform) printf("  [Uniform Start + HPWL Optimization]");
    else if (use_diffuse) printf("  [Diffusion Potential Spreading]");
    else if (use_spread) printf("  [Spread-First + CWTR Recovery]");
    else if (use_cwtr) printf("  [+CWTR Trust-Region Refinement]");
    if (use_bo) printf("  [Bayesian Lambda Opt]");
    printf("\n");
    printf("==========================================================\n\n");

    srand(42);

    // ---- Parse ----
    std::vector<Cell> cells;
    std::vector<Net> nets;
    float_t chip_xl, chip_yl, chip_xh, chip_yh;
    int num_terminals;

    if (!parse_ispd2005(bench_name, cells, nets, chip_xl, chip_yl, chip_xh, chip_yh, num_terminals)) {
        fprintf(stderr, "Error: Failed to parse benchmark\n");
        return 1;
    }

    PlacementResult result = {0, 0, 0, 0};
    PlacementResult cwtr_result = {0, 0, 0, 0};

    // ---- Apply initial placement strategy ----
    InitStrategy istrat = InitStrategy::EPLACE_QP;
    if (init_strategy == "uniform")      istrat = InitStrategy::UNIFORM_GRID;
    else if (init_strategy == "random")  istrat = InitStrategy::RANDOM;
    else if (init_strategy == "row")     istrat = InitStrategy::ROW_ORDERED;
    else if (init_strategy == "scaled")  istrat = InitStrategy::SCALED_EPLACE;
    else if (init_strategy == "gaussian") istrat = InitStrategy::GAUSSIAN_CENTER;
    else if (init_strategy == "quadrant") istrat = InitStrategy::QUADRANT_BINS;
    else if (init_strategy == "eplace")  istrat = InitStrategy::EPLACE_QP;

    if (istrat != InitStrategy::EPLACE_QP && !use_uniform) {
        // Override cell positions with chosen initial placement
        // (skip for --uniform because it has its own init in solve_uniform_start)
        printf("\n--- Applying initial placement: %s ---\n", init_strategy.c_str());
        std::string applied = apply_initial_placement(
            istrat, cells, nets, chip_xl, chip_yl, chip_xh, chip_yh, bench_name);
        printf("  Initial HPWL = %.1f\n", compute_hpwl_only(cells, nets));
    }

    const bool alternate_solver = use_spread || use_diffuse || use_uniform;
    if (alternate_solver && visualization.enabled) {
        visualization_initialize(visualization, cells, chip_xl, chip_yl, chip_xh, chip_yh);
        visualization_snapshot(visualization, 0, "initial", compute_hpwl_only(cells, nets),
                               0.0f, 0.0f, 0.0f, cells);
    }

    if (use_spread) {
        printf("══════════════════════════════════════════════════\n");
        printf("  SPREAD-FIRST: Aggressive Spread + CWTR Recovery\n");
        printf("══════════════════════════════════════════════════\n");
        result = solve_spread_first(cells, nets, chip_xl, chip_yl, chip_xh, chip_yh, 0.01f);
    } else if (use_diffuse) {
        printf("══════════════════════════════════════════════════\n");
        printf("  DIFFUSION: Smoothed Potential Field Spreading\n");
        printf("══════════════════════════════════════════════════\n");
        result = solve_diffusion_spread(cells, nets, chip_xl, chip_yl, chip_xh, chip_yh, 0.01f);
    } else if (use_uniform) {
        printf("══════════════════════════════════════════════════\n");
        printf("  UNIFORM: Uniform Start + HPWL Optimization\n");
        printf("══════════════════════════════════════════════════\n");
        result = solve_uniform_start(cells, nets, chip_xl, chip_yl, chip_xh, chip_yh);
    } else {
        // ---- Phase A: Subgradient Solver ----
        SolverConfig cfg;
        cfg.max_iters = max_iters;
        cfg.step_scale = 1.0f;
        cfg.step_min = 0.001f;
        cfg.step_max = 100.0f;
        cfg.momentum = 0.85f;
        cfg.lambda_init = 0.01f;
        cfg.lambda_growth = 1.15f;
        cfg.lambda_decay = 0.97f;
        cfg.target_overflow = 0.20f;
        cfg.grad_clip = 1000.0f;
        cfg.p0_iters = std::min(p0_iters, max_iters);
        cfg.coarse_iters = std::min(p1_iters, max_iters - cfg.p0_iters);
        cfg.medium_iters = std::min(p2_iters, max_iters - cfg.p0_iters - cfg.coarse_iters);
        cfg.p1_lambda_init = p1_lambda;
        cfg.use_bayesian_opt = use_bo;
        cfg.lambda_update_interval = lambda_interval;
        cfg.lambda_increase_factor = lambda_up;
        cfg.lambda_decrease_factor = lambda_down;
        cfg.lambda_ratio_hysteresis = lambda_hysteresis;
        cfg.lambda_control_max = lambda_max;
        cfg.overflow_baseline = overflow_baseline;
        cfg.overflow_limit = overflow_limit > 0.0f
            ? overflow_limit : overflow_baseline;
        cfg.target_density_override = target_density_override;
        cfg.use_adam = use_adam;
        cfg.use_lambda_guard = use_lambda_guard;
        cfg.use_trajectory_lambda = use_trajectory_lambda;
        cfg.lambda_guard_low = lambda_guard_low;
        cfg.lambda_guard_floor = lambda_guard_floor;
        cfg.lambda_guard_interval = lambda_guard_interval;
        cfg.lambda_trajectory_horizon = lambda_trajectory_horizon;
        cfg.lambda_trajectory_target = lambda_trajectory_target;
        cfg.fine_momentum = fine_momentum;
        cfg.fine_cooldown_floor = fine_cooldown_floor;
        cfg.use_projected_hpwl = use_projected_hpwl;
        cfg.use_adaptive_restart = use_adaptive_restart;
        cfg.projection_overflow_margin = projection_overflow_margin;
        cfg.use_hpwl_bundle = use_hpwl_bundle;
        cfg.hpwl_bundle_size = hpwl_bundle_size;
        cfg.hpwl_bundle_interval = hpwl_bundle_interval;
        cfg.hpwl_bundle_adaptive_cuts = hpwl_bundle_adaptive_cuts;
        cfg.hpwl_bundle_min_interval = hpwl_bundle_min_interval;
        cfg.hpwl_bundle_cosine_trigger = hpwl_bundle_cosine_trigger;
        cfg.hpwl_bundle_prox_scale = hpwl_bundle_prox_scale;
        cfg.hpwl_bundle_current_mix = hpwl_bundle_current_mix;
        cfg.hpwl_bundle_extrema_mix = hpwl_bundle_extrema_mix;
        cfg.hpwl_bundle_mix_gain = hpwl_bundle_mix_gain;
        cfg.hpwl_bundle_mix_decay_steps = hpwl_bundle_mix_decay_steps;
        cfg.visualization = visualization;

        printf("\n══════════════════════════════════════════════════\n");
        printf("  PHASE A: Subgradient Method with Polyak Step\n");
        if (use_short_edge) printf("  [Short-Edge Density Gradient Enabled]\n");
        printf("══════════════════════════════════════════════════\n");
        set_short_edge_only(use_short_edge);
        result = solve_placement(cells, nets, chip_xl, chip_yl, chip_xh, chip_yh, cfg);

        // ---- Phase B: CWTR Refinement (optional) ----
        if (use_cwtr) {
            printf("\n══════════════════════════════════════════════════\n");
            printf("  PHASE B: Component-wise Trust Region (CWTR)\n");
            printf("══════════════════════════════════════════════════\n");
            CWTRConfig cwtr_cfg;
            cwtr_cfg.max_sweeps = 1;
            cwtr_cfg.candidates_per_dim = 5;
            cwtr_cfg.trust_radius_initial = 2.0f;
            cwtr_cfg.overflow_relax = 1.05f;
            cwtr_cfg.hpwl_improve_min = 0.9999f;
            cwtr_cfg.swap_attempts = 0;
            cwtr_cfg.target_density = target_density_override > 0.0f
                ? target_density_override : 0.80f;
            cwtr_cfg.overflow_limit = 0.06f;
            cwtr_result = solve_cwtr(cells, nets, chip_xl, chip_yl, chip_xh, chip_yh, cwtr_cfg);
        }
    }

    if (alternate_solver && visualization.enabled) {
        visualization_snapshot(visualization, result.iterations, "final", result.hpwl,
                               result.overflow, 0.0f, 0.0f, cells);
    }

    // ---- Legalize ----
    std::string scl_file = bench_name + ".scl";
    legalize_placement(cells, chip_xl, chip_yl, chip_xh, chip_yh, scl_file);

    // ---- Write results ----
    std::string out_pl = bench_name + ".nsp.pl";
    write_placement(out_pl, cells, chip_xl, chip_yl, chip_xh, chip_yh);

    // ---- Summary ----
    printf("\n==========================================================\n");
    printf("  NSP Global Placement Complete\n");
    printf("==========================================================\n");
    printf("  Benchmark:        %s\n", bench_name.c_str());
    printf("  Movable cells:    %zu\n", cells.size() - num_terminals);
    printf("  Nets:             %zu\n", nets.size());

    if (use_spread || use_diffuse || use_uniform) {
        printf("  ───────────────────────────────────────────\n");
        printf("  %s Result:\n", use_diffuse ? "Diffusion" : "Spread-First");
        printf("    HPWL:           %.2f\n", result.hpwl);
        printf("    Overflow:       %.6f\n", result.overflow);
        printf("    Iter/Moves:     %d\n", result.iterations);
        printf("    Runtime:        %.2f sec\n", result.runtime_sec);
        printf("  ───────────────────────────────────────────\n");
        printf("  HPWL vs ePlace IP: %.1f%%\n",
               (44135944.0f - result.hpwl) / 44135944.0f * 100);
    } else {
        printf("  ───────────────────────────────────────────\n");
        printf("  Phase A (Subgrad):\n");
        printf("    HPWL:           %.2f\n", result.hpwl);
        printf("    Overflow:       %.6f\n", result.overflow);
        printf("    Iterations:     %d\n", result.iterations);
        printf("    Runtime:        %.2f sec\n", result.runtime_sec);

        if (use_cwtr) {
            printf("  ───────────────────────────────────────────\n");
            printf("  Phase B (CWTR):\n");
            printf("    HPWL:           %.2f\n", cwtr_result.hpwl);
            printf("    Overflow:       %.6f\n", cwtr_result.overflow);
            printf("    Moves:          %d\n", cwtr_result.iterations);
            printf("    Runtime:        %.2f sec\n", cwtr_result.runtime_sec);
            printf("  ───────────────────────────────────────────\n");
            printf("  Total HPWL improvement: %.1f%%\n",
                   (44135944.0f - cwtr_result.hpwl) / 44135944.0f * 100);
        }
    }
    printf("  Output:           %s\n", out_pl.c_str());

    return 0;
}
