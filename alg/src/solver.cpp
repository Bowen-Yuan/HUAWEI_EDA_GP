#include "types.h"
#include "bayesian_lambda.h"
#include "visualization.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <fstream>

float_t compute_hpwl_and_gradient(const std::vector<Cell>&, const std::vector<Net>&,
                                   std::vector<float_t>&, std::vector<float_t>&);
float_t compute_hpwl_only(const std::vector<Cell>&, const std::vector<Net>&);
void density_init(BinGrid&, float_t, float_t, float_t, float_t,
                  const std::vector<Cell>&, float_t);
float_t density_compute_and_gradient(const BinGrid&, const std::vector<Cell>&,
                                      std::vector<float_t>&, std::vector<float_t>&,
                                      std::vector<float_t>&);
float_t density_evaluate(const BinGrid&, const std::vector<Cell>&);

// ============================================================================
// Non-smooth Global Placement Solver
// Augmented-Lagrangian-style subgradient method:
//   Phase 1 (HPWL-only):    λ=0, pure WL optimization from IP
//   Phase 2 (WL+spread):    λ grows if overflow > target
//   Phase 3 (fine-tune):    balanced WL-density with decaying λ
// ============================================================================

PlacementResult solve_placement(
    std::vector<Cell>& cells,
    std::vector<Net>& nets,
    float_t chip_xl, float_t chip_yl,
    float_t chip_xh, float_t chip_yh,
    const SolverConfig& cfg)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    const int N = (int)cells.size();
    const int M = (int)nets.size();

    std::vector<int> movable_ids;
    movable_ids.reserve(N);
    for (int i = 0; i < N; ++i)
        if (!cells[i].is_terminal)
            movable_ids.push_back(i);
    const int Nm = (int)movable_ids.size();

    printf("[Solver] %d movable cells, %d nets\n", Nm, M);

    double total_cell_area = 0.0;
    for (int i : movable_ids) total_cell_area += (double)cells[i].area;
    double chip_area = (double)(chip_xh - chip_xl) * (double)(chip_yh - chip_yl);
    float_t target_density = (float_t)(total_cell_area / chip_area * 0.75);

    printf("[Solver] design_density=%.4f target=%.4f\n",
           total_cell_area / chip_area, target_density);

    std::vector<float_t> grad_x(N, 0.0f);
    std::vector<float_t> grad_y(N, 0.0f);
    std::vector<float_t> mom_x(N, 0.0f);
    std::vector<float_t> mom_y(N, 0.0f);

    // ---- Bin grids ----
    BinGrid grid_fine, grid_medium, grid_coarse;
    density_init(grid_fine, chip_xl, chip_yl, chip_xh, chip_yh, cells, target_density);
    visualization_initialize(cfg.visualization, cells, chip_xl, chip_yl, chip_xh, chip_yh);
    visualization_snapshot(cfg.visualization, 0, "initial", compute_hpwl_only(cells, nets),
                           0.0f, 0.0f, 0.0f, cells);

    auto make_coarse = [&](BinGrid& g, int factor) {
        g = grid_fine;
        g.nx = std::max(4, grid_fine.nx / factor);
        g.ny = std::max(4, grid_fine.ny / factor);
        g.bin_w = (chip_xh - chip_xl) / g.nx;
        g.bin_h = (chip_yh - chip_yl) / g.ny;
        g.density.assign(g.nx * g.ny, 0.0f);
        g.overflow.assign(g.nx * g.ny, 0.0f);
    };
    make_coarse(grid_medium, 2);
    make_coarse(grid_coarse, 4);

    printf("[Solver] Grids: C=%dx%d M=%dx%d F=%dx%d\n",
           grid_coarse.nx, grid_coarse.ny,
           grid_medium.nx, grid_medium.ny,
           grid_fine.nx, grid_fine.ny);

    int max_bins = std::max({(int)grid_coarse.density.size(),
                             (int)grid_medium.density.size(),
                             (int)grid_fine.density.size()});
    std::vector<float_t> bin_buf(max_bins, 0.0f);

    float_t chip_w = chip_xh - chip_xl, chip_h = chip_yh - chip_yl;
    float_t hpwl = 0.0f, dpen = 0.0f;
    float_t hpwl_best = 1e30f;
    int total_iters = 0;
    const float_t RMS_EPS = 1e-6f;

    // ---- CSV convergence log ----
    // Logs: iter, phase, hpwl, objective, avg_overflow_real, dpen, lambda, g_rms, lr
    std::ofstream csv_log("nsp_convergence.csv");
    csv_log << "iter,phase,hpwl,objective,avg_overflow,dpen,lambda,g_rms,lr\n";

    // Helper: compute real avg_overflow from bin density data
    auto compute_real_overflow = [](const BinGrid& grid, const std::vector<float_t>& bin_area_sum,
                                     int nbins) -> float_t {
        float_t bin_area = grid.bin_w * grid.bin_h;
        float_t sum_over = 0.0f;
        int count = 0;
        for (int b = 0; b < nbins; ++b) {
            float_t rho = bin_area_sum[b] / bin_area;
            float_t over = rho - grid.target_density;
            if (over > 0.0f) {
                sum_over += over;
                ++count;
            }
        }
        return count > 0 ? sum_over / nbins : 0.0f;
    };

    // ================================================================
    // Phase 0: HPWL-only optimization (preserve QP quality)
    // ================================================================
    int n0 = cfg.p0_iters;
    printf("\n--- Phase 0: HPWL-only (%d iters, λ=0) ---\n", n0);
    float_t lambda = 0.0f;
    float_t base_step = (chip_w + chip_h) * 0.005f;

    for (int iter = 0; iter < n0; ++iter, ++total_iters) {
        hpwl = compute_hpwl_and_gradient(cells, nets, grad_x, grad_y);

        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            double gn = (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
            gsum2 += gn;
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        float_t progress = (float_t)iter / n0;
        float_t lr = base_step / std::sqrt(progress * 2.0f + 1.0f);

        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;
            gx = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gx));
            gy = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gy));

            float_t beta = 0.85f;
            mom_x[i] = beta * mom_x[i] + (1.0f - beta) * gx;
            mom_y[i] = beta * mom_y[i] + (1.0f - beta) * gy;

            cells[i].x -= lr * mom_x[i];
            cells[i].y -= lr * mom_y[i];
            float_t m = cells[i].width;
            cells[i].x = std::max(chip_xl - m, std::min(chip_xh + m, cells[i].x));
            cells[i].y = std::max(chip_yl - m, std::min(chip_yh + m, cells[i].y));
        }

        if (hpwl < hpwl_best) hpwl_best = hpwl;

        csv_log << iter << ",0," << hpwl << "," << hpwl << ",0,0,0,"
                << g_rms << "," << lr << "\n";
        visualization_snapshot(cfg.visualization, total_iters + 1, "hpwl", hpwl,
                               0.0f, 0.0f, 0.0f, cells);

        if (iter % 30 == 0 || iter == n0 - 1)
            printf("  [p0%3d] HPWL=%.1f lr=%.2f g_rms=%.3f\n",
                   iter, hpwl, lr, g_rms);
    }

    // ================================================================
    // Phase 1: Gentle density introduction with adaptive λ
    // ================================================================
    int n1 = cfg.coarse_iters;
    lambda = cfg.p1_lambda_init;
    base_step = (chip_w + chip_h) * 0.008f;
    float_t target_overflow = 0.10f;

    // Initialize Bayesian Optimizer for λ (if enabled)
    LambdaBayesianOpt lambda_bo;
    if (cfg.use_bayesian_opt) {
        LambdaBOConfig bo_cfg;
        bo_cfg.length_scale  = 1.5;
        bo_cfg.signal_std    = 0.3;
        bo_cfg.noise_std     = 0.03;
        bo_cfg.penalty_weight = 5.0;
        bo_cfg.window_size   = 35;
        bo_cfg.refit_interval = 5;
        bo_cfg.lambda_min    = 1e-4;
        bo_cfg.lambda_max    = 20.0;
        bo_cfg.grid_points   = 200;
        bo_cfg.ei_xi         = 0.02;
        lambda_bo = LambdaBayesianOpt(bo_cfg);
        printf("[BO] Bayesian optimization enabled for λ selection\n");
        printf("     GP(Matern5/2) + EI, refit every %d iters, window=%d\n",
               bo_cfg.refit_interval, bo_cfg.window_size);
    }

    printf("\n--- Phase 1: Coarse Spread (%d iters) ---\n", n1);
    for (int iter = 0; iter < n1; ++iter, ++total_iters) {
        const BinGrid& grid = grid_coarse;
        bin_buf.resize(grid.nx * grid.ny);

        hpwl = compute_hpwl_and_gradient(cells, nets, grad_x, grad_y);
        dpen = density_compute_and_gradient(grid, cells, grad_x, grad_y, bin_buf);

        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            double gn = (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
            gsum2 += gn;
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        float_t progress = (float_t)iter / n1;
        float_t lr = base_step / std::sqrt(progress * 4.0f + 1.0f);

        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;
            gx = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gx));
            gy = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gy));

            float_t beta = cfg.momentum;
            mom_x[i] = beta * mom_x[i] + (1.0f - beta) * gx;
            mom_y[i] = beta * mom_y[i] + (1.0f - beta) * gy;

            cells[i].x -= lr * mom_x[i];
            cells[i].y -= lr * mom_y[i];
            float_t m = cells[i].width;
            cells[i].x = std::max(chip_xl - m, std::min(chip_xh + m, cells[i].x));
            cells[i].y = std::max(chip_yl - m, std::min(chip_yh + m, cells[i].y));
        }

        if (hpwl < hpwl_best) hpwl_best = hpwl;

        float_t avg_over = (dpen > 0) ? std::sqrt(dpen / (grid.nx * grid.ny)) : 0.0f;

        float_t real_overflow = compute_real_overflow(grid, bin_buf, grid.nx * grid.ny);

        if (cfg.use_bayesian_opt) {
            lambda = (float_t)lambda_bo.step((double)lambda, (double)hpwl,
                                              (double)real_overflow, (double)target_overflow);
        } else {
            if (avg_over > target_overflow * 1.5f)
                lambda = std::min(20.0f, lambda * 1.08f);
            else if (avg_over > target_overflow)
                lambda = std::min(20.0f, lambda * 1.03f);
            else
                lambda = std::max(0.0005f, lambda * 0.98f);
        }

        float_t objective = hpwl + lambda * dpen;

        csv_log << iter << ",1," << hpwl << "," << objective << ","
                << real_overflow << "," << dpen << "," << lambda << ","
                << g_rms << "," << lr << "\n";
        visualization_snapshot(cfg.visualization, total_iters + 1, "coarse", hpwl,
                               real_overflow, dpen, lambda, cells);

        if (iter % 30 == 0 || iter == n1 - 1)
            printf("  [c%3d] HPWL=%.1f Dpen=%.1f λ=%.4f lr=%.2f g_rms=%.3f ov=%.4f\n",
                   iter, hpwl, dpen, lambda, lr, g_rms, avg_over);
    }

    // Phase transition: force BO to re-explore
    if (cfg.use_bayesian_opt) {
        float_t real_ov = compute_real_overflow(grid_coarse, bin_buf, grid_coarse.nx * grid_coarse.ny);
        lambda = (float_t)lambda_bo.force_select((double)lambda, (double)hpwl,
                                                  (double)real_ov, (double)target_overflow);
        lambda_bo.reset();
        printf("  [BO] Phase 1→2 transition: λ reset to %.4f\n", lambda);
    }

    // ================================================================
    // Phase 2: Medium grid balance
    // ================================================================
    int n2 = cfg.medium_iters;
    base_step = (chip_w + chip_h) * 0.004f;

    printf("\n--- Phase 2: Medium Balance (%d iters) ---\n", n2);
    for (int iter = 0; iter < n2; ++iter, ++total_iters) {
        const BinGrid& grid = grid_medium;
        bin_buf.resize(grid.nx * grid.ny);

        hpwl = compute_hpwl_and_gradient(cells, nets, grad_x, grad_y);
        dpen = density_compute_and_gradient(grid, cells, grad_x, grad_y, bin_buf);

        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            double gn = (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
            gsum2 += gn;
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        float_t progress = (float_t)iter / n2;
        float_t lr = base_step / std::sqrt(progress * 4.0f + 1.0f);

        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;
            gx = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gx));
            gy = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gy));

            float_t beta = 0.80f;
            mom_x[i] = beta * mom_x[i] + (1.0f - beta) * gx;
            mom_y[i] = beta * mom_y[i] + (1.0f - beta) * gy;

            cells[i].x -= lr * mom_x[i];
            cells[i].y -= lr * mom_y[i];
            float_t m = cells[i].width;
            cells[i].x = std::max(chip_xl - m, std::min(chip_xh + m, cells[i].x));
            cells[i].y = std::max(chip_yl - m, std::min(chip_yh + m, cells[i].y));
        }

        if (hpwl < hpwl_best) hpwl_best = hpwl;

        float_t avg_over = (dpen > 0) ? std::sqrt(dpen / (grid.nx * grid.ny)) : 0.0f;

        float_t real_overflow = compute_real_overflow(grid, bin_buf, grid.nx * grid.ny);

        if (cfg.use_bayesian_opt) {
            lambda = (float_t)lambda_bo.step((double)lambda, (double)hpwl,
                                              (double)real_overflow, (double)target_overflow);
        } else {
            if (avg_over > target_overflow * 1.5f)
                lambda = std::min(20.0f, lambda * 1.05f);
            else if (avg_over > target_overflow)
                lambda = std::min(20.0f, lambda * 1.01f);
            else
                lambda = std::max(0.0005f, lambda * 0.995f);
        }

        float_t objective = hpwl + lambda * dpen;

        csv_log << iter << ",2," << hpwl << "," << objective << ","
                << real_overflow << "," << dpen << "," << lambda << ","
                << g_rms << "," << lr << "\n";
        visualization_snapshot(cfg.visualization, total_iters + 1, "medium", hpwl,
                               real_overflow, dpen, lambda, cells);

        if (iter % 50 == 0 || iter == n2 - 1)
            printf("  [m%3d] HPWL=%.1f Dpen=%.1f λ=%.4f lr=%.2f g_rms=%.3f ov=%.4f\n",
                   iter, hpwl, dpen, lambda, lr, g_rms, avg_over);
    }

    // Phase transition: force BO to re-explore
    if (cfg.use_bayesian_opt) {
        lambda = (float_t)lambda_bo.force_select((double)lambda, (double)hpwl,
                                                  (double)0.3, (double)target_overflow);
        lambda_bo.reset();
        printf("  [BO] Phase 2→3 transition: λ reset to %.4f\n", lambda);
    }

    // ================================================================
    // Phase 3: Fine grid refinement
    // ================================================================
    int n3 = cfg.max_iters - total_iters;
    base_step = (chip_w + chip_h) * 0.002f;

    printf("\n--- Phase 3: Fine Refinement (%d iters) ---\n", n3);
    for (int iter = 0; iter < n3; ++iter, ++total_iters) {
        const BinGrid& grid = grid_fine;
        bin_buf.resize(grid.nx * grid.ny);

        hpwl = compute_hpwl_and_gradient(cells, nets, grad_x, grad_y);
        dpen = density_compute_and_gradient(grid, cells, grad_x, grad_y, bin_buf);

        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            double gn = (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
            gsum2 += gn;
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        float_t progress = (float_t)iter / n3;
        float_t lr = base_step / std::sqrt(progress * 8.0f + 1.0f);

        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;
            gx = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gx));
            gy = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gy));

            float_t beta = 0.70f;
            mom_x[i] = beta * mom_x[i] + (1.0f - beta) * gx;
            mom_y[i] = beta * mom_y[i] + (1.0f - beta) * gy;

            cells[i].x -= lr * mom_x[i];
            cells[i].y -= lr * mom_y[i];
            cells[i].x = std::max(chip_xl, std::min(chip_xh, cells[i].x));
            cells[i].y = std::max(chip_yl, std::min(chip_yh, cells[i].y));
        }

        if (hpwl < hpwl_best) hpwl_best = hpwl;

        float_t avg_over = (dpen > 0) ? std::sqrt(dpen / (grid.nx * grid.ny)) : 0.0f;

        float_t real_overflow = compute_real_overflow(grid, bin_buf, grid.nx * grid.ny);

        if (cfg.use_bayesian_opt) {
            lambda = (float_t)lambda_bo.step((double)lambda, (double)hpwl,
                                              (double)real_overflow, (double)target_overflow);
        } else {
            if (avg_over > target_overflow)
                lambda = std::min(5.0f, lambda * 1.01f);
            else
                lambda = std::max(0.0001f, lambda * 0.998f);
        }

        float_t objective = hpwl + lambda * dpen;

        csv_log << iter << ",3," << hpwl << "," << objective << ","
                << real_overflow << "," << dpen << "," << lambda << ","
                << g_rms << "," << lr << "\n";
        visualization_snapshot(cfg.visualization, total_iters + 1, "fine", hpwl,
                               real_overflow, dpen, lambda, cells);

        if (iter % 100 == 0 || iter == n3 - 1)
            printf("  [f%3d] HPWL=%.1f Dpen=%.1f λ=%.5f lr=%.2f g_rms=%.3f ov=%.4f\n",
                   iter, hpwl, dpen, lambda, lr, g_rms, avg_over);
    }

    csv_log.close();
    printf("[Solver] Convergence log saved to nsp_convergence.csv\n");

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    float_t final_hpwl = compute_hpwl_only(cells, nets);
    float_t final_overflow = density_evaluate(grid_fine, cells);

    printf("\n==========================================\n");
    printf("  Done: %d iters, %.2f sec\n", total_iters, elapsed);
    printf("  HPWL final=%.1f  best=%.1f  overflow=%.4f\n",
           final_hpwl, hpwl_best, final_overflow);
    printf("==========================================\n");

    return { final_hpwl, final_overflow, total_iters, elapsed };
}
