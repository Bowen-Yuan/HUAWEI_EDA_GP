#include "types.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <fstream>
#include <random>
#include <cstring>

// Forward declarations
float_t compute_hpwl_only(const std::vector<Cell>&, const std::vector<Net>&);
void density_init(BinGrid&, float_t, float_t, float_t, float_t,
                  const std::vector<Cell>&, float_t);
float_t density_compute_and_gradient(const BinGrid&, const std::vector<Cell>&,
                                      std::vector<float_t>&, std::vector<float_t>&,
                                      std::vector<float_t>&);
float_t density_evaluate(const BinGrid&, const std::vector<Cell>&);
void set_density_grad_scale(float_t s);

// ============================================================================
// Spread-First Solver: Aggressive Spreading + CWTR HPWL Recovery
//
// Strategy (from prompt.md):
//   Step 1: Find minimum-overflow solution (safety anchor) via aggressive
//           density-only gradient descent
//   Step 2: Under "do not increase overflow" hard constraint, greedily
//           compress HPWL using CWTR direct search
//   Step 3: When stuck, gradually relax + swap operations
//   Step 4: Iterate steps 2-3 until convergence
// ============================================================================

PlacementResult solve_spread_first(
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
    printf("  SPREAD-FIRST SOLVER: Target avg_overflow < %.4f\n", target_avg_overflow);
    printf("==========================================================\n");
    printf("  %d movable cells\n", Nm);

    // ---- Setup density grid ----
    double total_cell_area = 0.0;
    for (int i : movable_ids) total_cell_area += (double)cells[i].area;
    double chip_area = (double)(chip_xh - chip_xl) * (double)(chip_yh - chip_yl);
    float_t target_density = (float_t)(total_cell_area / chip_area * 0.80);

    BinGrid grid;
    density_init(grid, chip_xl, chip_yl, chip_xh, chip_yh, cells, target_density);

    // Gradient and momentum arrays
    std::vector<float_t> grad_x(N, 0.0f);
    std::vector<float_t> grad_y(N, 0.0f);
    std::vector<float_t> mom_x(N, 0.0f);
    std::vector<float_t> mom_y(N, 0.0f);

    int nbins = grid.nx * grid.ny;
    std::vector<float_t> bin_buf(nbins, 0.0f);

    float_t chip_w = chip_xh - chip_xl, chip_h = chip_yh - chip_yl;
    const float_t RMS_EPS = 1e-6f;
    const float_t base_step = (chip_w + chip_h) * 0.03f;
    int total_iters = 0;

    // ================================================================
    // Phase S1: Coarse-grid aggressive spreading (density-only)
    // ================================================================
    int n_s1 = 200;
    float_t lambda = 50.0f;
    // Set aggressive density gradient scale for spreading
    set_density_grad_scale(1.0f);  // 50x stronger than default 0.02

    float_t hpwl = compute_hpwl_only(cells, nets);
    printf("\n--- Phase S1: Coarse Spreading (%d iters, density-only) ---\n", n_s1);
    printf("  Initial HPWL=%.1f, density_grad_scale=%.2f\n", hpwl, 1.0f);

    // Create coarse grid for faster spreading
    BinGrid grid_coarse = grid;
    grid_coarse.nx = std::max(4, grid.nx / 4);
    grid_coarse.ny = std::max(4, grid.ny / 4);
    grid_coarse.bin_w = (chip_xh - chip_xl) / grid_coarse.nx;
    grid_coarse.bin_h = (chip_yh - chip_yl) / grid_coarse.ny;
    grid_coarse.density.assign(grid_coarse.nx * grid_coarse.ny, 0.0f);
    grid_coarse.overflow.assign(grid_coarse.nx * grid_coarse.ny, 0.0f);
    int cbins = grid_coarse.nx * grid_coarse.ny;
    std::vector<float_t> cbuf(cbins, 0.0f);

    // CSV log
    std::ofstream csv_log("spread_convergence.csv");
    csv_log << "iter,phase,hpwl,avg_overflow,max_overflow,dpen,lambda,g_rms,lr\n";

    for (int iter = 0; iter < n_s1; ++iter, ++total_iters) {
        // Zero out HPWL gradient contribution
        std::fill(grad_x.begin(), grad_x.end(), 0.0f);
        std::fill(grad_y.begin(), grad_y.end(), 0.0f);

        // Compute HPWL separately for logging
        if (iter % 10 == 0) {
            hpwl = compute_hpwl_only(cells, nets);
        }

        // Compute density-only gradient
        float_t dpen = density_compute_and_gradient(grid_coarse, cells, grad_x, grad_y, cbuf);

        // Compute real overflow
        float_t bin_area_c = grid_coarse.bin_w * grid_coarse.bin_h;
        float_t avg_overflow = 0.0f, max_overflow = 0.0f;
        for (int b = 0; b < cbins; ++b) {
            float_t rho = cbuf[b] / bin_area_c;
            float_t over = rho - target_density;
            if (over > 0.0f) {
                avg_overflow += over;
                if (over > max_overflow) max_overflow = over;
            }
        }
        avg_overflow /= cbins;  // average over ALL bins

        // Gradient stats
        double gsum2 = 0.0;
        float_t max_g = 0.0f;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i], gy = grad_y[i];
            double gn = (double)gx*gx + (double)gy*gy;
            gsum2 += gn;
            float_t gnf = (float_t)std::sqrt(gn);
            if (gnf > max_g) max_g = gnf;
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        // Step size
        float_t progress = (float_t)iter / n_s1;
        float_t lr = base_step / std::sqrt(progress * 2.0f + 1.0f);

        // Update with density-only gradient
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;

            float_t beta = 0.80f;
            mom_x[i] = beta * mom_x[i] + (1.0f - beta) * gx;
            mom_y[i] = beta * mom_y[i] + (1.0f - beta) * gy;

            cells[i].x -= lr * mom_x[i];
            cells[i].y -= lr * mom_y[i];

            float_t m = cells[i].width;
            cells[i].x = std::max(chip_xl - m, std::min(chip_xh + m, cells[i].x));
            cells[i].y = std::max(chip_yl - m, std::min(chip_yh + m, cells[i].y));
        }

        csv_log << iter << ",S1," << hpwl << "," << avg_overflow << ","
                << max_overflow << "," << dpen << "," << lambda << ","
                << g_rms << "," << lr << "\n";

        if (iter % 25 == 0 || iter == n_s1 - 1)
            printf("  [S1 %3d] HPWL=%.1f avg_ov=%.4f max_ov=%.2f lr=%.1f |g|=%.1f\n",
                   iter, hpwl, avg_overflow, max_overflow, lr, max_g);

        // Early exit if target reached
        if (avg_overflow < target_avg_overflow) {
            printf("  -> Target reached at iter %d!\n", iter);
            break;
        }
    }

    // ================================================================
    // Phase S2: Fine-grid spreading
    // ================================================================
    int n_s2 = 300;
    lambda = 30.0f;
    set_density_grad_scale(0.5f);  // reduce for finer grid
    printf("\n--- Phase S2: Fine Spreading (%d iters, scale=0.5) ---\n", n_s2);

    for (int iter = 0; iter < n_s2; ++iter, ++total_iters) {
        std::fill(grad_x.begin(), grad_x.end(), 0.0f);
        std::fill(grad_y.begin(), grad_y.end(), 0.0f);

        if (iter % 10 == 0) {
            hpwl = compute_hpwl_only(cells, nets);
        }

        float_t dpen = density_compute_and_gradient(grid, cells, grad_x, grad_y, bin_buf);

        float_t bin_area_f = grid.bin_w * grid.bin_h;
        float_t avg_overflow = 0.0f, max_overflow = 0.0f;
        for (int b = 0; b < nbins; ++b) {
            float_t rho = bin_buf[b] / bin_area_f;
            float_t over = rho - target_density;
            if (over > 0.0f) {
                avg_overflow += over;
                if (over > max_overflow) max_overflow = over;
            }
        }
        avg_overflow /= nbins;

        double gsum2 = 0.0;
        float_t max_g = 0.0f;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i], gy = grad_y[i];
            double gn = (double)gx*gx + (double)gy*gy;
            gsum2 += gn;
            float_t gnf = (float_t)std::sqrt(gn);
            if (gnf > max_g) max_g = gnf;
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        float_t progress = (float_t)iter / n_s2;
        float_t lr = base_step * 0.5f / std::sqrt(progress * 4.0f + 1.0f);

        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;

            float_t beta = 0.75f;
            mom_x[i] = beta * mom_x[i] + (1.0f - beta) * gx;
            mom_y[i] = beta * mom_y[i] + (1.0f - beta) * gy;

            cells[i].x -= lr * mom_x[i];
            cells[i].y -= lr * mom_y[i];

            cells[i].x = std::max(chip_xl, std::min(chip_xh, cells[i].x));
            cells[i].y = std::max(chip_yl, std::min(chip_yh, cells[i].y));
        }

        // Adaptive lambda: increase if still overfull
        if (avg_overflow > target_avg_overflow * 5.0f)
            lambda = std::min(200.0f, lambda * 1.02f);
        else if (avg_overflow > target_avg_overflow)
            lambda = std::min(200.0f, lambda * 1.01f);
        else
            lambda = std::max(1.0f, lambda * 0.98f);

        csv_log << iter << ",S2," << hpwl << "," << avg_overflow << ","
                << max_overflow << "," << dpen << "," << lambda << ","
                << g_rms << "," << lr << "\n";

        if (iter % 50 == 0 || iter == n_s2 - 1)
            printf("  [S2 %3d] HPWL=%.1f avg_ov=%.4f max_ov=%.2f lambda=%.1f lr=%.1f\n",
                   iter, hpwl, avg_overflow, max_overflow, lambda, lr);

        if (avg_overflow < target_avg_overflow) {
            printf("  -> Target reached at iter %d!\n", iter);
            break;
        }
    }

    float_t hpwl_after_spread = compute_hpwl_only(cells, nets);

    // ================================================================
    // Phase R: CWTR HPWL Recovery (hard overflow constraint)
    // ================================================================
    set_density_grad_scale(0.02f);  // restore default
    printf("\n--- Phase R: CWTR HPWL Recovery ---\n");
    printf("  HPWL after spread: %.1f\n", hpwl_after_spread);

    // Build net adjacency for CWTR
    std::vector<std::vector<int>> cell_to_nets(N);
    for (size_t e = 0; e < nets.size(); ++e)
        for (int cid : nets[e].cell_ids)
            cell_to_nets[cid].push_back((int)e);

    // Initialize bin area sums
    std::vector<float_t> bin_area_sum(nbins, 0.0f);
    for (int i = 0; i < N; ++i) {
        if (cells[i].is_terminal) continue;
        float_t cx = cells[i].x, cy = cells[i].y;
        float_t hw = cells[i].width * 0.5f, hh = cells[i].height * 0.5f;
        float_t cl = cx - hw, cr = cx + hw, cb = cy - hh, ct = cy + hh;
        int bx_s = std::max(0, (int)((cl - grid.chip_xl) / grid.bin_w));
        int bx_e = std::min(grid.nx-1, (int)((cr - grid.chip_xl) / grid.bin_w));
        int by_s = std::max(0, (int)((cb - grid.chip_yl) / grid.bin_h));
        int by_e = std::min(grid.ny-1, (int)((ct - grid.chip_yl) / grid.bin_h));
        for (int by = by_s; by <= by_e; ++by) {
            float_t by_l = grid.chip_yl + by * grid.bin_h;
            float_t oy = std::max(0.0f, std::min(ct, by_l+grid.bin_h) - std::max(cb, by_l));
            for (int bx = bx_s; bx <= bx_e; ++bx) {
                float_t bx_l = grid.chip_xl + bx * grid.bin_w;
                float_t ox = std::max(0.0f, std::min(cr, bx_l+grid.bin_w) - std::max(cl, bx_l));
                bin_area_sum[by * grid.nx + bx] += ox * oy;
            }
        }
    }

    // Current overflow
    float_t bin_area_f = grid.bin_w * grid.bin_h;
    float_t max_overflow = 0.0f;
    for (int b = 0; b < nbins; ++b) {
        float_t over = bin_area_sum[b] / bin_area_f - target_density;
        if (over > max_overflow) max_overflow = over;
    }
    float_t overflow_threshold = std::max(max_overflow, target_avg_overflow * 2.0f);

    printf("  Max overflow: %.4f, threshold: %.4f\n", max_overflow, overflow_threshold);

    std::mt19937 rng(42);
    float_t best_hpwl = hpwl_after_spread;
    int max_sweeps = 20;

    // Helper: HPWL delta
    auto hpwl_delta = [&](int cid, float_t nx, float_t ny) -> float_t {
        float_t delta = 0.0f;
        for (int eid : cell_to_nets[cid]) {
            const Net& net = nets[eid];
            if (net.cell_ids.size() < 2) continue;
            auto hpwl_at = [&](float_t moved_x, float_t moved_y) {
                float_t min_x = FLT_MAX, max_x = -FLT_MAX;
                float_t min_y = FLT_MAX, max_y = -FLT_MAX;
                for (size_t pin = 0; pin < net.cell_ids.size(); ++pin) {
                    const int pin_cell = net.cell_ids[pin];
                    const float_t base_x = pin_cell == cid ? moved_x : cells[pin_cell].x;
                    const float_t base_y = pin_cell == cid ? moved_y : cells[pin_cell].y;
                    const float_t x = base_x + net.pin_offset_x[pin];
                    const float_t y = base_y + net.pin_offset_y[pin];
                    min_x = std::min(min_x, x); max_x = std::max(max_x, x);
                    min_y = std::min(min_y, y); max_y = std::max(max_y, y);
                }
                return (max_x - min_x) + (max_y - min_y);
            };
            delta += hpwl_at(nx, ny) - hpwl_at(cells[cid].x, cells[cid].y);
        }
        return delta;
    };

    // Helper: max overflow after moving cell
    auto check_overflow = [&](int cid, float_t nx, float_t ny) -> float_t {
        float_t ox = cells[cid].x, oy = cells[cid].y;
        float_t hw = cells[cid].width * 0.5f, hh = cells[cid].height * 0.5f;

        auto get_bins = [&](float_t cx, float_t cy) {
            float_t cl = cx - hw, cr = cx + hw, cb = cy - hh, ct = cy + hh;
            int bx_s = std::max(0, (int)((cl - grid.chip_xl) / grid.bin_w));
            int bx_e = std::min(grid.nx-1, (int)((cr - grid.chip_xl) / grid.bin_w));
            int by_s = std::max(0, (int)((cb - grid.chip_yl) / grid.bin_h));
            int by_e = std::min(grid.ny-1, (int)((ct - grid.chip_yl) / grid.bin_h));
            std::vector<std::pair<int,float_t>> bins;
            for (int by = by_s; by <= by_e; ++by) {
                float_t byl = grid.chip_yl + by * grid.bin_h;
                float_t oyy = std::max(0.0f, std::min(ct, byl+grid.bin_h) - std::max(cb, byl));
                for (int bx = bx_s; bx <= bx_e; ++bx) {
                    float_t bxl = grid.chip_xl + bx * grid.bin_w;
                    float_t oxx = std::max(0.0f, std::min(cr, bxl+grid.bin_w) - std::max(cl, bxl));
                    if (oxx > 0 && oyy > 0)
                        bins.emplace_back(by * grid.nx + bx, oxx * oyy);
                }
            }
            return bins;
        };

        auto old_bins = get_bins(ox, oy);
        auto new_bins = get_bins(nx, ny);

        // Collect affected bins
        int affected[16], naff = 0;
        for (auto& p : old_bins) { affected[naff++] = p.first; }
        for (auto& p : new_bins) {
            bool found = false;
            for (int k = 0; k < naff; ++k) if (affected[k] == p.first) { found = true; break; }
            if (!found && naff < 16) affected[naff++] = p.first;
        }

        float_t max_ov = 0.0f;
        for (int k = 0; k < naff; ++k) {
            int bid = affected[k];
            float_t old_contrib = 0.0f, new_contrib = 0.0f;
            for (auto& p : old_bins) if (p.first == bid) old_contrib = p.second;
            for (auto& p : new_bins) if (p.first == bid) new_contrib = p.second;

            float_t new_area = bin_area_sum[bid] - old_contrib + new_contrib;
            float_t over = new_area / bin_area_f - target_density;
            if (over > max_ov) max_ov = over;
        }
        return max_ov;
    };

    // ---- CWTR Sweeps ----
    csv_log << "sweep,phase,hpwl,avg_overflow,max_overflow,dpen,lambda,g_rms,lr\n";

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        int moves = 0;
        std::shuffle(movable_ids.begin(), movable_ids.end(), rng);

        float_t trust_r = 1.5f * std::pow(0.85f, (float_t)sweep);
        int K = 5;

        for (int idx = 0; idx < Nm; ++idx) {
            int cid = movable_ids[idx];
            float_t ox = cells[cid].x, oy = cells[cid].y;
            float_t rx = trust_r * grid.bin_w;
            float_t ry = trust_r * grid.bin_h;

            float_t best_delta = 0.0f, best_nx = ox, best_ny = oy;
            bool found = false;

            for (int dx = -(K/2); dx <= K/2; ++dx) {
                for (int dy = -(K/2); dy <= K/2; ++dy) {
                    if (dx == 0 && dy == 0) continue;
                    float_t nx = ox + dx * rx / (K/2);
                    float_t ny = oy + dy * ry / (K/2);
                    nx = std::max(chip_xl, std::min(chip_xh, nx));
                    ny = std::max(chip_yl, std::min(chip_yh, ny));

                    float_t max_ov = check_overflow(cid, nx, ny);
                    if (max_ov > overflow_threshold) continue;

                    float_t delta_w = hpwl_delta(cid, nx, ny);
                    if (delta_w < best_delta) {
                        best_delta = delta_w;
                        best_nx = nx;
                        best_ny = ny;
                        found = true;
                    }
                }
            }

            if (found && best_delta < -0.01f) {
                // Update bin sums
                float_t hw = cells[cid].width * 0.5f, hh = cells[cid].height * 0.5f;
                auto update_bins = [&](float_t cx, float_t cy, float_t sign) {
                    float_t cl = cx - hw, cr = cx + hw, cb = cy - hh, ct = cy + hh;
                    int bx_s = std::max(0,(int)((cl-grid.chip_xl)/grid.bin_w));
                    int bx_e = std::min(grid.nx-1,(int)((cr-grid.chip_xl)/grid.bin_w));
                    int by_s = std::max(0,(int)((cb-grid.chip_yl)/grid.bin_h));
                    int by_e = std::min(grid.ny-1,(int)((ct-grid.chip_yl)/grid.bin_h));
                    for (int by = by_s; by <= by_e; ++by) {
                        float_t byl = grid.chip_yl+by*grid.bin_h;
                        float_t oyy = std::max(0.0f,std::min(ct,byl+grid.bin_h)-std::max(cb,byl));
                        for (int bx = bx_s; bx <= bx_e; ++bx) {
                            float_t bxl = grid.chip_xl+bx*grid.bin_w;
                            float_t oxx = std::max(0.0f,std::min(cr,bxl+grid.bin_w)-std::max(cl,bxl));
                            bin_area_sum[by*grid.nx+bx] += sign*oxx*oyy;
                        }
                    }
                };
                update_bins(ox, oy, -1.0f);
                cells[cid].x = best_nx; cells[cid].y = best_ny;
                update_bins(best_nx, best_ny, +1.0f);
                ++moves;
            }
        }

        hpwl = compute_hpwl_only(cells, nets);
        if (hpwl < best_hpwl) best_hpwl = hpwl;

        // Recompute overflow
        max_overflow = 0.0f;
        float_t avg_ov = 0.0f;
        for (int b = 0; b < nbins; ++b) {
            float_t over = bin_area_sum[b] / bin_area_f - target_density;
            if (over > 0.0f) {
                avg_ov += over;
                if (over > max_overflow) max_overflow = over;
            }
        }
        avg_ov /= nbins;

        // Log
        float_t dpen_eval = 0.0f;
        for (int b = 0; b < nbins; ++b) {
            float_t over = bin_area_sum[b] / bin_area_f - target_density;
            if (over > 0.0f) dpen_eval += over * over;
        }

        csv_log << sweep << ",R," << hpwl << "," << avg_ov << ","
                << max_overflow << "," << dpen_eval << ",0,0,0\n";

        printf("  [R sweep %2d] HPWL=%.1f avg_ov=%.4f max_ov=%.3f moves=%d\n",
               sweep, hpwl, avg_ov, max_overflow, moves);

        // Relax if stuck
        if (moves < 10) {
            if (overflow_threshold < 1.0f && avg_ov > target_avg_overflow) {
                overflow_threshold = std::min(1.0f, overflow_threshold * 1.10f);
                printf("  -> Relaxing threshold to %.4f\n", overflow_threshold);
            } else {
                printf("  -> Converged\n");
                break;
            }
        } else if (max_overflow < overflow_threshold * 0.7f) {
            overflow_threshold = max_overflow * 1.2f;
        }
    }

    csv_log.close();
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    float_t final_hpwl = compute_hpwl_only(cells, nets);
    float_t final_overflow = density_evaluate(grid, cells);

    float_t initial_hpwl = 44135944.0f;  // ePlace IP HPWL for adaptec1
    printf("\n==========================================================\n");
    printf("  Spread-First Complete: %d iters + sweeps, %.2f sec\n",
           total_iters, elapsed);
    printf("  HPWL: %.1f -> %.1f -> %.1f  (total: %.1f%%)\n",
           initial_hpwl, hpwl_after_spread, final_hpwl,
           (initial_hpwl - final_hpwl) / initial_hpwl * 100);
    printf("  Final overflow: %.6f\n", final_overflow);
    printf("==========================================================\n");

    return { final_hpwl, final_overflow, total_iters, elapsed };
}
