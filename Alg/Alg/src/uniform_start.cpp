#include "types.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <fstream>
#include <random>

// Forward declarations
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
// Uniform-Start Solver: start from uniform placement (~0% overflow),
// then optimize HPWL while monitoring overflow.
// ============================================================================

PlacementResult solve_uniform_start(
    std::vector<Cell>& cells,
    std::vector<Net>& nets,
    float_t chip_xl, float_t chip_yl,
    float_t chip_xh, float_t chip_yh)
{
    auto t_start = std::chrono::high_resolution_clock::now();
    const int N = (int)cells.size();

    std::vector<int> movable_ids;
    movable_ids.reserve(N);
    for (int i = 0; i < N; ++i)
        if (!cells[i].is_terminal)
            movable_ids.push_back(i);
    const int Nm = (int)movable_ids.size();

    printf("\n==========================================================\n");
    printf("  UNIFORM-START SOLVER\n");
    printf("==========================================================\n");

    float_t chip_w = chip_xh - chip_xl, chip_h = chip_yh - chip_yl;

    // ---- Initialize cells uniformly across chip area ----
    // Sort cells by width for better packing
    std::sort(movable_ids.begin(), movable_ids.end(), [&](int a, int b) {
        return cells[a].width > cells[b].width;
    });

    // Place cells in a grid pattern
    int cols = (int)std::sqrt((double)Nm * chip_w / chip_h);
    cols = std::max(1, cols);
    float_t cell_spacing_x = chip_w / cols;
    float_t cell_spacing_y = chip_h / ((Nm + cols - 1) / cols);

    for (int idx = 0; idx < Nm; ++idx) {
        int i = movable_ids[idx];
        int row = idx / cols;
        int col = idx % cols;
        cells[i].x = chip_xl + (col + 0.5f) * cell_spacing_x;
        cells[i].y = chip_yl + (row + 0.5f) * cell_spacing_y;

        // Clamp
        cells[i].x = std::max(chip_xl + cells[i].width*0.5f,
                              std::min(chip_xh - cells[i].width*0.5f, cells[i].x));
        cells[i].y = std::max(chip_yl + cells[i].height*0.5f,
                              std::min(chip_yh - cells[i].height*0.5f, cells[i].y));
    }

    float_t hpwl = compute_hpwl_only(cells, nets);
    printf("  Uniform-start HPWL=%.1f\n", hpwl);

    // Setup density grid
    double total_cell_area = 0.0;
    for (int i : movable_ids) total_cell_area += (double)cells[i].area;
    double chip_area = (double)(chip_xh - chip_xl) * (double)(chip_yh - chip_yl);
    float_t target_density = (float_t)(total_cell_area / chip_area * 0.75);

    BinGrid grid;
    density_init(grid, chip_xl, chip_yl, chip_xh, chip_yh, cells, target_density);

    // Gradient arrays
    std::vector<float_t> grad_x(N, 0.0f);
    std::vector<float_t> grad_y(N, 0.0f);
    std::vector<float_t> mom_x(N, 0.0f);
    std::vector<float_t> mom_y(N, 0.0f);

    int nbins = grid.nx * grid.ny;
    std::vector<float_t> bin_buf(nbins, 0.0f);

    const float_t RMS_EPS = 1e-6f;
    const float_t base_step = (chip_w + chip_h) * 0.008f;
    int total_iters = 0;

    // CSV log
    std::ofstream csv_log("uniform_convergence.csv");
    csv_log << "iter,hpwl,avg_overflow,max_overflow,dpen,lambda,g_rms,lr\n";

    // ---- Phase 1: HPWL optimization with very weak density ----
    float_t lambda = 0.001f;
    int n1 = 300;
    float_t hpwl_best = 1e30f;

    printf("\n--- Phase U1: HPWL Optimization (lambda=%.4f, %d iters) ---\n", lambda, n1);
    for (int iter = 0; iter < n1; ++iter, ++total_iters) {
        hpwl = compute_hpwl_and_gradient(cells, nets, grad_x, grad_y);
        float_t dpen = density_compute_and_gradient(grid, cells, grad_x, grad_y, bin_buf);

        // Overflow stats
        float_t bin_area = grid.bin_w * grid.bin_h;
        float_t avg_ov = 0.0f, max_ov = 0.0f;
        for (int b = 0; b < nbins; ++b) {
            float_t over = bin_buf[b] / bin_area - target_density;
            if (over > 0.0f) {
                avg_ov += over;
                if (over > max_ov) max_ov = over;
            }
        }
        avg_ov /= nbins;

        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            gsum2 += (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        float_t progress = (float_t)iter / n1;
        float_t lr = base_step / std::sqrt(progress * 3.0f + 1.0f);

        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;

            float_t beta = 0.85f;
            mom_x[i] = beta * mom_x[i] + (1.0f-beta) * gx;
            mom_y[i] = beta * mom_y[i] + (1.0f-beta) * gy;

            cells[i].x -= lr * mom_x[i];
            cells[i].y -= lr * mom_y[i];
            cells[i].x = std::max(chip_xl, std::min(chip_xh, cells[i].x));
            cells[i].y = std::max(chip_yl, std::min(chip_yh, cells[i].y));
        }

        if (hpwl < hpwl_best) hpwl_best = hpwl;

        // Adaptive lambda: increase if overflow grows
        if (avg_ov > 0.15f)
            lambda = std::min(0.5f, lambda * 1.05f);
        else if (avg_ov > 0.05f)
            lambda = std::min(0.5f, lambda * 1.01f);

        csv_log << iter << "," << hpwl << "," << avg_ov << "," << max_ov
                << "," << dpen << "," << lambda << "," << g_rms << "," << lr << "\n";

        if (iter % 50 == 0 || iter == n1-1)
            printf("  [U1 %3d] HPWL=%.1f avg_ov=%.4f max_ov=%.2f lambda=%.4f lr=%.1f\n",
                   iter, hpwl, avg_ov, max_ov, lambda, lr);
    }

    // ---- Phase 2: More HPWL optimization with moderate density control ----
    int n2 = 300;
    lambda = std::min(0.5f, lambda);
    printf("\n--- Phase U2: Balanced (%d iters, lambda=%.4f) ---\n", n2, lambda);

    for (int iter = 0; iter < n2; ++iter, ++total_iters) {
        hpwl = compute_hpwl_and_gradient(cells, nets, grad_x, grad_y);
        float_t dpen = density_compute_and_gradient(grid, cells, grad_x, grad_y, bin_buf);

        float_t bin_area = grid.bin_w * grid.bin_h;
        float_t avg_ov = 0.0f, max_ov = 0.0f;
        for (int b = 0; b < nbins; ++b) {
            float_t over = bin_buf[b] / bin_area - target_density;
            if (over > 0.0f) { avg_ov += over; if (over > max_ov) max_ov = over; }
        }
        avg_ov /= nbins;

        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            gsum2 += (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        float_t progress = (float_t)iter / n2;
        float_t lr = base_step * 0.5f / std::sqrt(progress * 5.0f + 1.0f);

        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;
            float_t beta = 0.80f;
            mom_x[i] = beta * mom_x[i] + (1.0f-beta) * gx;
            mom_y[i] = beta * mom_y[i] + (1.0f-beta) * gy;
            cells[i].x -= lr * mom_x[i];
            cells[i].y -= lr * mom_y[i];
            cells[i].x = std::max(chip_xl, std::min(chip_xh, cells[i].x));
            cells[i].y = std::max(chip_yl, std::min(chip_yh, cells[i].y));
        }

        if (hpwl < hpwl_best) hpwl_best = hpwl;

        if (avg_ov > 0.10f) lambda = std::min(1.0f, lambda * 1.02f);
        else lambda = std::max(0.0001f, lambda * 0.995f);

        csv_log << iter << "," << hpwl << "," << avg_ov << "," << max_ov
                << "," << dpen << "," << lambda << "," << g_rms << "," << lr << "\n";

        if (iter % 50 == 0 || iter == n2-1)
            printf("  [U2 %3d] HPWL=%.1f avg_ov=%.4f max_ov=%.2f lambda=%.4f lr=%.1f\n",
                   iter, hpwl, avg_ov, max_ov, lambda, lr);
    }

    csv_log.close();
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    float_t final_hpwl = compute_hpwl_only(cells, nets);
    float_t final_overflow = density_evaluate(grid, cells);

    printf("\n==========================================================\n");
    printf("  Uniform-Start Complete: %d iters, %.2f sec\n", total_iters, elapsed);
    printf("  HPWL: best=%.1f, final=%.1f, overflow=%.4f\n",
           hpwl_best, final_hpwl, final_overflow);
    printf("==========================================================\n");

    return { final_hpwl, final_overflow, total_iters, elapsed };
}
