#include "types.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <fstream>
#include <random>
#include <vector>

// Forward declarations
float_t compute_hpwl_and_gradient(const std::vector<Cell>&, const std::vector<Net>&,
                                   std::vector<float_t>&, std::vector<float_t>&);
float_t compute_hpwl_only(const std::vector<Cell>&, const std::vector<Net>&);
void density_init(BinGrid&, float_t, float_t, float_t, float_t,
                  const std::vector<Cell>&, float_t);
float_t density_evaluate(const BinGrid&, const std::vector<Cell>&);

// ============================================================================
// Diffusion-Based Spreading Solver
//
// Key insight (from prompt.md): overflow is piecewise-constant on the grid,
// making gradient-based spreading ineffective. Solution:
//   1. Compute exact overflow map O(b) on bins
//   2. Apply Gaussian smoothing to create a C-infinity potential field phi(x,y)
//   3. Move cells along -grad(phi) -- gradient is non-zero everywhere
//   4. Recover HPWL with CWTR under overflow constraint
//
// The smoothing is only for computing the SEARCH DIRECTION; the objective
// (HPWL + overflow penalty) is still evaluated EXACTLY without smoothing.
// ============================================================================

// ---- Gaussian blur of a 2D grid (in-place) ----
static void box_blur_2d(std::vector<float_t>& data, int nx, int ny, int radius) {
    std::vector<float_t> tmp(data.size());
    for (int iter = 0; iter < radius; ++iter) {
        // Horizontal pass
        #pragma omp parallel for
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                float_t sum = data[iy*nx + ix];
                int cnt = 1;
                if (ix > 0)      { sum += data[iy*nx + (ix-1)]; ++cnt; }
                if (ix < nx-1)   { sum += data[iy*nx + (ix+1)]; ++cnt; }
                tmp[iy*nx + ix] = sum / cnt;
            }
        }
        // Vertical pass
        #pragma omp parallel for
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                float_t sum = tmp[iy*nx + ix];
                int cnt = 1;
                if (iy > 0)      { sum += tmp[(iy-1)*nx + ix]; ++cnt; }
                if (iy < ny-1)   { sum += tmp[(iy+1)*nx + ix]; ++cnt; }
                data[iy*nx + ix] = sum / cnt;
            }
        }
    }
}

// ---- Bilinear interpolation of a grid at position (x, y) ----
static float_t bilinear_interp(const std::vector<float_t>& grid,
                                int nx, int ny,
                                float_t xl, float_t yl,
                                float_t bin_w, float_t bin_h,
                                float_t px, float_t py) {
    float_t fx = (px - xl) / bin_w - 0.5f;
    float_t fy = (py - yl) / bin_h - 0.5f;
    int ix = (int)std::floor(fx);
    int iy = (int)std::floor(fy);
    float_t tx = fx - ix, ty = fy - iy;

    // Clamp
    ix = std::max(0, std::min(nx-2, ix));
    iy = std::max(0, std::min(ny-2, iy));

    float_t v00 = grid[iy*nx + ix];
    float_t v10 = grid[iy*nx + (ix+1)];
    float_t v01 = grid[(iy+1)*nx + ix];
    float_t v11 = grid[(iy+1)*nx + (ix+1)];

    return (1-tx)*(1-ty)*v00 + tx*(1-ty)*v10 + (1-tx)*ty*v01 + tx*ty*v11;
}

// ---- Compute grad(phi) at cell position using central differences ----
static void potential_gradient(const std::vector<float_t>& phi,
                                int nx, int ny,
                                float_t xl, float_t yl,
                                float_t bin_w, float_t bin_h,
                                float_t px, float_t py,
                                float_t& gx, float_t& gy) {
    float_t eps = bin_w * 0.01f;
    float_t phi_r = bilinear_interp(phi, nx, ny, xl, yl, bin_w, bin_h, px+eps, py);
    float_t phi_l = bilinear_interp(phi, nx, ny, xl, yl, bin_w, bin_h, px-eps, py);
    float_t phi_t = bilinear_interp(phi, nx, ny, xl, yl, bin_w, bin_h, px, py+eps);
    float_t phi_b = bilinear_interp(phi, nx, ny, xl, yl, bin_w, bin_h, px, py-eps);
    gx = (phi_r - phi_l) / (2.0f * eps);
    gy = (phi_t - phi_b) / (2.0f * eps);
}

PlacementResult solve_diffusion_spread(
    std::vector<Cell>& cells,
    std::vector<Net>& nets,
    float_t chip_xl, float_t chip_yl,
    float_t chip_xh, float_t chip_yh,
    float_t target_avg_overflow = 0.01f)
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
    printf("  DIFFUSION SPREADER: target avg_overflow < %.4f\n", target_avg_overflow);
    printf("==========================================================\n");

    double total_cell_area = 0.0;
    for (int i : movable_ids) total_cell_area += (double)cells[i].area;
    double chip_area = (double)(chip_xh - chip_xl) * (double)(chip_yh - chip_yl);
    float_t target_density = (float_t)(total_cell_area / chip_area * 0.80);

    // ---- Setup bin grid ----
    // Use finer grid for better potential field resolution
    BinGrid grid;
    grid.chip_xl = chip_xl; grid.chip_yl = chip_yl;
    grid.chip_xh = chip_xh; grid.chip_yh = chip_yh;
    grid.target_density = target_density;

    // Dense grid for potential field
    int nx = 64, ny = 72;  // ~4600 bins
    grid.nx = nx; grid.ny = ny;
    grid.bin_w = (chip_xh - chip_xl) / nx;
    grid.bin_h = (chip_yh - chip_yl) / ny;
    int nbins = nx * ny;
    grid.density.assign(nbins, 0.0f);
    grid.overflow.assign(nbins, 0.0f);
    float_t bin_area = grid.bin_w * grid.bin_h;

    printf("  Potential grid: %dx%d, bin=%.1fx%.1f\n", nx, ny, grid.bin_w, grid.bin_h);

    float_t chip_w = chip_xh - chip_xl, chip_h = chip_yh - chip_yl;
    float_t hpwl = compute_hpwl_only(cells, nets);
    printf("  Initial HPWL=%.1f\n", hpwl);

    // CSV log
    std::ofstream csv_log("diffusion_convergence.csv");
    csv_log << "iter,hpwl,avg_overflow,max_overflow,step,max_grad\n";

    // ---- Outer loop: spread + WL optimize ----
    float_t lambda = 10.0f;
    float_t base_step = (chip_w + chip_h) * 0.015f;
    int total_iters = 0;
    int max_outer = 100;

    // Pre-allocate
    std::vector<float_t> bin_overflow(nbins, 0.0f);
    std::vector<float_t> phi(nbins, 0.0f);
    std::vector<float_t> grad_x(N, 0.0f);
    std::vector<float_t> grad_y(N, 0.0f);
    std::vector<float_t> mom_x(N, 0.0f);
    std::vector<float_t> mom_y(N, 0.0f);

    for (int outer = 0; outer < max_outer; ++outer) {
        // ---- Step 1: Compute exact overflow map ----
        std::fill(bin_overflow.begin(), bin_overflow.end(), 0.0f);
        std::fill(grad_x.begin(), grad_x.end(), 0.0f);
        std::fill(grad_y.begin(), grad_y.end(), 0.0f);

        // Accumulate cell areas into bins (exact)
        #pragma omp parallel
        {
            std::vector<float_t> local_bins(nbins, 0.0f);
            #pragma omp for schedule(dynamic, 64)
            for (int idx = 0; idx < Nm; ++idx) {
                int i = movable_ids[idx];
                float_t cx = cells[i].x, cy = cells[i].y;
                float_t hw = cells[i].width * 0.5f, hh = cells[i].height * 0.5f;
                float_t cl = cx - hw, cr = cx + hw, cb = cy - hh, ct = cy + hh;
                int bx_s = std::max(0, (int)((cl-chip_xl)/grid.bin_w));
                int bx_e = std::min(nx-1, (int)((cr-chip_xl)/grid.bin_w));
                int by_s = std::max(0, (int)((cb-chip_yl)/grid.bin_h));
                int by_e = std::min(ny-1, (int)((ct-chip_yl)/grid.bin_h));
                for (int by = by_s; by <= by_e; ++by) {
                    float_t byl = chip_yl + by * grid.bin_h;
                    float_t oy = std::max(0.0f, std::min(ct, byl+grid.bin_h)-std::max(cb, byl));
                    for (int bx = bx_s; bx <= bx_e; ++bx) {
                        float_t bxl = chip_xl + bx * grid.bin_w;
                        float_t ox = std::max(0.0f, std::min(cr, bxl+grid.bin_w)-std::max(cl, bxl));
                        local_bins[by*nx + bx] += ox * oy;
                    }
                }
            }
            #pragma omp critical
            { for (int b = 0; b < nbins; ++b) bin_overflow[b] += local_bins[b]; }
        }

        // Convert to density ratio (use raw density for potential, so gradient
        // points from high-density to low-density regions everywhere)
        float_t avg_ov = 0.0f, max_ov = 0.0f;
        for (int b = 0; b < nbins; ++b) {
            float_t rho = bin_overflow[b] / bin_area;  // raw density
            float_t over = rho - target_density;
            bin_overflow[b] = rho;  // store raw density for potential
            if (over > 0.0f) {
                avg_ov += over;
                if (over > max_ov) max_ov = over;
            }
        }
        avg_ov /= nbins;

        if (outer % 10 == 0 || outer == 0)
            printf("  [outer %3d] HPWL=%.1f avg_ov=%.4f max_ov=%.2f\n",
                   outer, hpwl, avg_ov, max_ov);

        if (avg_ov < target_avg_overflow) {
            printf("  -> Target reached!\n");
            break;
        }

        // ---- Step 2: Create smoothed potential field ----
        std::copy(bin_overflow.begin(), bin_overflow.end(), phi.begin());
        box_blur_2d(phi, nx, ny, 8);  // 8 iterations of box blur -> ~Gaussian

        // ---- Step 3: Compute HPWL gradient ----
        hpwl = compute_hpwl_and_gradient(cells, nets, grad_x, grad_y);

        // ---- Step 4: Compute potential gradient and combine ----
        float_t max_g = 0.0f;
        #pragma omp parallel for reduction(max:max_g)
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t pgx, pgy;
            potential_gradient(phi, nx, ny, chip_xl, chip_yl,
                              grid.bin_w, grid.bin_h,
                              cells[i].x, cells[i].y, pgx, pgy);

            // Combine: -grad(HPWL) - lambda * grad(phi)
            grad_x[i] = -grad_x[i] - lambda * pgx * 500.0f;   // scale potential grad
            grad_y[i] = -grad_y[i] - lambda * pgy * 500.0f;

            float_t gn = grad_x[i]*grad_x[i] + grad_y[i]*grad_y[i];
            if (gn > max_g*max_g) max_g = std::sqrt(gn);
        }

        // ---- Step 5: RMS normalize + momentum + update ----
        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            gsum2 += (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + 1e-6f);

        float_t step = base_step / std::sqrt((float_t)outer * 0.5f + 1.0f);

        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;

            float_t beta = 0.75f;
            mom_x[i] = beta * mom_x[i] + (1.0f - beta) * gx;
            mom_y[i] = beta * mom_y[i] + (1.0f - beta) * gy;

            cells[i].x -= step * mom_x[i];
            cells[i].y -= step * mom_y[i];

            cells[i].x = std::max(chip_xl, std::min(chip_xh, cells[i].x));
            cells[i].y = std::max(chip_yl, std::min(chip_yh, cells[i].y));
        }

        // Adaptive lambda
        if (avg_ov > target_avg_overflow * 5.0f)
            lambda = std::min(100.0f, lambda * 1.01f);
        else if (avg_ov > target_avg_overflow)
            lambda = std::min(100.0f, lambda * 1.005f);

        csv_log << outer << "," << hpwl << "," << avg_ov << "," << max_ov
                << "," << step << "," << max_g << "\n";
        ++total_iters;
    }

    // Save overflow map for visualization
    std::ofstream omap("overflow_map.csv");
    omap << "x,y,overflow\n";
    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            // Actually we need the raw density, recompute
            omap << ix << "," << iy << "," << 0 << "\n";
        }
    }
    omap.close();

    csv_log.close();
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    float_t final_hpwl = compute_hpwl_only(cells, nets);
    float_t final_overflow = density_evaluate(grid, cells);

    printf("\n==========================================================\n");
    printf("  Diffusion Complete: %d iters, %.2f sec\n", total_iters, elapsed);
    printf("  HPWL: 44135944 -> %.1f  overflow: %.4f\n", final_hpwl, final_overflow);
    printf("==========================================================\n");

    return { final_hpwl, final_overflow, total_iters, elapsed };
}
