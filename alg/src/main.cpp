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
        printf("  --short-edge  : Zero out less-effective density gradient direction\n");
        printf("  --visualize [DIR] : export optional snapshots (default visualizations)\n");
        printf("  --visualize-every N : snapshot interval (default 10)\n");
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
    int max_iters = 800;         // Total NSP budget
    float_t p1_lambda = 0.005f;  // Phase 1 initial λ
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
            cwtr_cfg.max_sweeps = 15;
            cwtr_cfg.candidates_per_dim = 5;
            cwtr_cfg.trust_radius_initial = 2.0f;
            cwtr_cfg.overflow_relax = 1.05f;
            cwtr_cfg.hpwl_improve_min = 0.9999f;
            cwtr_cfg.swap_attempts = 500;
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
