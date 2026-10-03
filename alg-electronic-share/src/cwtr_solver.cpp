#include "types.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <fstream>
#include <random>
#include <set>

// ============================================================================
// Component-wise Trust Region (CWTR) Solver
//
// Based on the observation that overflow is piecewise-constant on the bin grid:
//   - Moving a cell within a bin: overflow unchanged (flat region)
//   - Crossing a bin boundary: overflow jumps (discontinuity)
//
// The subgradient method struggles with this because gradients vanish in flat
// regions and oscillate at boundaries. CWTR addresses this by:
//
//   1. TRUST REGION: each cell i can move freely within a local neighborhood
//      of bins, as long as overflow does not exceed a threshold.
//   2. GREEDY HPWL: within the trust region, find the position that minimizes
//      HPWL contribution of cell i (only nets connected to i).
//   3. CONSTRAINT RELAXATION: when no cell can improve HPWL without violating
//      the overflow constraint, gradually relax the threshold.
//   4. SWAP OPERATIONS: swap positions of two cells to restructure configuration.
//   5. ITERATE: repeat 2-4 until convergence.
// ============================================================================

// Forward declarations
float_t compute_hpwl_only(const std::vector<Cell>&, const std::vector<Net>&);
void density_init(BinGrid&, float_t, float_t, float_t, float_t,
                  const std::vector<Cell>&, float_t);
float_t density_evaluate(const BinGrid&, const std::vector<Cell>&);
void density_build_map(const BinGrid&, const std::vector<Cell>&,
                       std::vector<float_t>&);

// ---------------------------------------------------------------------------
// Incremental HPWL computation for a single cell
// ---------------------------------------------------------------------------
struct NetAdj {
    std::vector<std::vector<int>> cell_to_nets;  // cell_id → net_ids
};

static NetAdj build_adjacency(const std::vector<Cell>& cells,
                               const std::vector<Net>& nets) {
    NetAdj adj;
    adj.cell_to_nets.resize(cells.size());
    for (size_t e = 0; e < nets.size(); ++e) {
        for (int cid : nets[e].cell_ids) {
            adj.cell_to_nets[cid].push_back((int)e);
        }
    }
    return adj;
}

// Precomputed bounds for a net without a specific cell
struct NetBounds {
    float_t min_x, max_x, min_y, max_y;
    int n_min_x, n_max_x;  // count of pins at extremum (for tie-breaking)
    int n_min_y, n_max_y;
};

// Compute bounds of a net WITHOUT cell `exclude_cid`
static NetBounds net_bounds_without(const Net& net, const std::vector<Cell>& cells,
                                     int exclude_cid) {
    NetBounds b;
    b.min_x = 1e30f; b.max_x = -1e30f;
    b.min_y = 1e30f; b.max_y = -1e30f;
    b.n_min_x = b.n_max_x = b.n_min_y = b.n_max_y = 0;

    for (int pid : net.cell_ids) {
        if (pid == exclude_cid) continue;
        float_t cx = cells[pid].x, cy = cells[pid].y;
        if (cx < b.min_x)      { b.min_x = cx; b.n_min_x = 1; }
        else if (cx == b.min_x) { b.n_min_x++; }
        if (cx > b.max_x)      { b.max_x = cx; b.n_max_x = 1; }
        else if (cx == b.max_x) { b.n_max_x++; }
        if (cy < b.min_y)      { b.min_y = cy; b.n_min_y = 1; }
        else if (cy == b.min_y) { b.n_min_y++; }
        if (cy > b.max_y)      { b.max_y = cy; b.n_max_y = 1; }
        else if (cy == b.max_y) { b.n_max_y++; }
    }
    return b;
}

// HPWL delta when moving cell `cid` to (new_x, new_y) — O(num_nets) per call
static float_t compute_hpwl_delta(int cid, float_t new_x, float_t new_y,
                                   const std::vector<Cell>& cells,
                                   const std::vector<Net>& nets,
                                   const NetAdj& adj) {
    float_t old_x = cells[cid].x, old_y = cells[cid].y;
    float_t delta = 0.0f;

    for (int eid : adj.cell_to_nets[cid]) {
        const Net& net = nets[eid];
        if (net.cell_ids.size() < 2) continue;

        // Compute bounds without this cell (once)
        NetBounds b = net_bounds_without(net, cells, cid);

        // If the net has no other pins, skip
        if (b.n_min_x == 0) continue;

        // HPWL before (with cell at old position)
        float_t old_min_x = std::min(b.min_x, old_x);
        float_t old_max_x = std::max(b.max_x, old_x);
        float_t old_min_y = std::min(b.min_y, old_y);
        float_t old_max_y = std::max(b.max_y, old_y);
        float_t hpwl_old = (old_max_x - old_min_x) + (old_max_y - old_min_y);

        // HPWL after (with cell at new position)
        float_t new_min_x = std::min(b.min_x, new_x);
        float_t new_max_x = std::max(b.max_x, new_x);
        float_t new_min_y = std::min(b.min_y, new_y);
        float_t new_max_y = std::max(b.max_y, new_y);
        float_t hpwl_new = (new_max_x - new_min_x) + (new_max_y - new_min_y);

        delta += hpwl_new - hpwl_old;
    }

    return delta;
}

// ---------------------------------------------------------------------------
// Incremental overflow computation for a single cell move
// ---------------------------------------------------------------------------
static float_t compute_overflow_after_move(
        int cid, float_t new_x, float_t new_y,
        const std::vector<Cell>& cells, const BinGrid& grid,
        const std::vector<float_t>& bin_area_sum,
        float_t current_excess_area) {
    float_t old_x = cells[cid].x, old_y = cells[cid].y;
    float_t hw = cells[cid].width * 0.5f, hh = cells[cid].height * 0.5f;

    auto bin_range = [&](float_t cx, float_t cy) {
        float_t cl = cx - hw, cr = cx + hw;
        float_t cb = cy - hh, ct = cy + hh;
        int bx_s = std::max(0, (int)((cl - grid.chip_xl) / grid.bin_w));
        int bx_e = std::min(grid.nx-1, (int)((cr - grid.chip_xl) / grid.bin_w));
        int by_s = std::max(0, (int)((cb - grid.chip_yl) / grid.bin_h));
        int by_e = std::min(grid.ny-1, (int)((ct - grid.chip_yl) / grid.bin_h));
        return std::make_tuple(bx_s, bx_e, by_s, by_e);
    };

    auto cell_area_in_bins = [&](float_t cx, float_t cy,
                                  int bx_s, int bx_e, int by_s, int by_e) {
        float_t cl = cx - hw, cr = cx + hw;
        float_t cb = cy - hh, ct = cy + hh;
        std::vector<std::pair<int,float_t>> bins;
        for (int by = by_s; by <= by_e; ++by) {
            float_t bin_yl = grid.chip_yl + by * grid.bin_h;
            float_t bin_yh = bin_yl + grid.bin_h;
            float_t oy = std::max(0.0f, std::min(ct, bin_yh) - std::max(cb, bin_yl));
            for (int bx = bx_s; bx <= bx_e; ++bx) {
                float_t bin_xl = grid.chip_xl + bx * grid.bin_w;
                float_t bin_xh = bin_xl + grid.bin_w;
                float_t ox = std::max(0.0f, std::min(cr, bin_xh) - std::max(cl, bin_xl));
                if (ox > 0 && oy > 0)
                    bins.emplace_back(by * grid.nx + bx, ox * oy);
            }
        }
        return bins;
    };

    auto [ox_s, ox_e, oy_s, oy_e] = bin_range(old_x, old_y);
    auto old_bins = cell_area_in_bins(old_x, old_y, ox_s, ox_e, oy_s, oy_e);

    auto [nx_s, nx_e, ny_s, ny_e] = bin_range(new_x, new_y);
    auto new_bins = cell_area_in_bins(new_x, new_y, nx_s, nx_e, ny_s, ny_e);

    // Collect all affected bins
    std::set<int> affected;
    for (auto& p : old_bins) affected.insert(p.first);
    for (auto& p : new_bins) affected.insert(p.first);

    float_t bin_area = grid.bin_w * grid.bin_h;
    float_t excess_area = current_excess_area;

    for (int bid : affected) {
        // Find old and new contributions
        float_t old_contrib = 0.0f, new_contrib = 0.0f;
        for (auto& p : old_bins) if (p.first == bid) old_contrib = p.second;
        for (auto& p : new_bins) if (p.first == bid) new_contrib = p.second;

        const float_t new_area = bin_area_sum[bid] - old_contrib + new_contrib;
        const float_t capacity = grid.target_density * bin_area;
        excess_area += std::max(0.0f, new_area - capacity) -
                       std::max(0.0f, bin_area_sum[bid] - capacity);
    }

    return std::max(0.0f, excess_area);
}

// ---------------------------------------------------------------------------
// Main CWTR optimization pass
// ---------------------------------------------------------------------------

PlacementResult solve_cwtr(
    std::vector<Cell>& cells,
    std::vector<Net>& nets,
    float_t chip_xl, float_t chip_yl,
    float_t chip_xh, float_t chip_yh,
    const CWTRConfig& cwtr_cfg = CWTRConfig{})
{
    auto t_start = std::chrono::high_resolution_clock::now();
    const std::vector<Cell> initial_cells = cells;

    const int N = (int)cells.size();
    std::vector<int> movable_ids;
    for (int i = 0; i < N; ++i)
        if (!cells[i].is_terminal)
            movable_ids.push_back(i);
    const int Nm = (int)movable_ids.size();

    printf("[CWTR] %d movable cells\n", Nm);

    // Build net adjacency
    printf("[CWTR] Building net adjacency...\n");
    NetAdj adj = build_adjacency(cells, nets);

    // Setup density grid
    const float_t target_density = cwtr_cfg.target_density;

    BinGrid grid;
    density_init(grid, chip_xl, chip_yl, chip_xh, chip_yh, cells, target_density);

    // Initialize bin area sums
    int nbins = grid.nx * grid.ny;
    std::vector<float_t> bin_area_sum;
    density_build_map(grid, cells, bin_area_sum);

    // Initial overflow
    float_t bin_area = grid.bin_w * grid.bin_h;
    float_t current_excess_area = 0.0f;
    for (int b = 0; b < nbins; ++b) {
        current_excess_area += std::max(0.0f,
            bin_area_sum[b] - target_density * bin_area);
    }
    float_t overflow_ratio = current_excess_area / grid.available_area;
    float_t overflow_threshold = cwtr_cfg.overflow_limit;
    float_t max_overflow = overflow_ratio;

    printf("[CWTR] target_density=%.4f initial_overflow=%.6f limit=%.6f\n",
           target_density, overflow_ratio, overflow_threshold);

    // CSV logging
    std::ofstream csv_log("cwtr_convergence.csv");
    csv_log << "sweep,hpwl,overflow,limit,moves_accepted\n";

    float_t initial_hpwl = compute_hpwl_only(cells, nets);
    float_t best_hpwl = initial_hpwl;
    printf("[CWTR] Initial HPWL=%.1f\n", initial_hpwl);

    std::mt19937 rng(42);
    int total_moves = 0;

    // ---- CWTR Main Loop ----
    for (int sweep = 0; sweep < cwtr_cfg.max_sweeps; ++sweep) {
        int moves_this_sweep = 0;

        // Shuffle cell order for fairness
        std::shuffle(movable_ids.begin(), movable_ids.end(), rng);

        float_t trust_radius = cwtr_cfg.trust_radius_initial
                               * std::pow(0.85f, (float_t)sweep);

        for (int idx = 0; idx < Nm; ++idx) {
            int cid = movable_ids[idx];
            float_t old_x = cells[cid].x, old_y = cells[cid].y;

            // Generate candidate positions in trust region
            float_t rx = trust_radius * grid.bin_w;
            float_t ry = trust_radius * grid.bin_h;
            int K = cwtr_cfg.candidates_per_dim;
            float_t best_hpwl_delta = 0.0f;
            float_t best_nx = old_x, best_ny = old_y;
            float_t best_excess_area = current_excess_area;

            for (int dx = -(K/2); dx <= K/2; ++dx) {
                for (int dy = -(K/2); dy <= K/2; ++dy) {
                    if (dx == 0 && dy == 0) continue;

                    float_t nx = old_x + dx * rx / (K/2);
                    float_t ny = old_y + dy * ry / (K/2);

                    // Clamp to chip area
                    nx = std::max(chip_xl, std::min(chip_xh, nx));
                    ny = std::max(chip_yl, std::min(chip_yh, ny));

                    // Check overflow constraint
                    float_t candidate_excess = compute_overflow_after_move(
                        cid, nx, ny, cells, grid, bin_area_sum,
                        current_excess_area);

                    if (candidate_excess / grid.available_area >
                        overflow_threshold + 1.0e-7f) continue;

                    // Compute HPWL delta
                    float_t hpwl_delta = compute_hpwl_delta(
                        cid, nx, ny, cells, nets, adj);

                    if (hpwl_delta < best_hpwl_delta) {
                        best_hpwl_delta = hpwl_delta;
                        best_nx = nx;
                        best_ny = ny;
                        best_excess_area = candidate_excess;
                    }
                }
            }

            // Accept move if HPWL improves
            if (best_hpwl_delta < -0.01f) {
                // Update bin area sums
                float_t hw = cells[cid].width * 0.5f;
                float_t hh = cells[cid].height * 0.5f;

                auto update_bins = [&](float_t cx, float_t cy, float_t sign) {
                    float_t cl = cx - hw, cr = cx + hw;
                    float_t cb = cy - hh, ct = cy + hh;
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
                            bin_area_sum[by*grid.nx + bx] += sign * ox * oy;
                        }
                    }
                };

                update_bins(old_x, old_y, -1.0f);  // remove from old
                cells[cid].x = best_nx;
                cells[cid].y = best_ny;
                update_bins(best_nx, best_ny, +1.0f);  // add to new
                current_excess_area = best_excess_area;

                ++moves_this_sweep;
                ++total_moves;
            }
        }

        // Evaluate
        float_t cur_hpwl = compute_hpwl_only(cells, nets);

        // Recompute the exact baseline overflow to eliminate incremental drift.
        current_excess_area = 0.0f;
        for (int b = 0; b < nbins; ++b) {
            current_excess_area += std::max(0.0f,
                bin_area_sum[b] - target_density * bin_area);
        }
        overflow_ratio = current_excess_area / grid.available_area;

        csv_log << sweep << "," << cur_hpwl << "," << overflow_ratio << ","
                << overflow_threshold << "," << moves_this_sweep << "\n";

        printf("  [sweep %2d] HPWL=%.1f overflow=%.6f limit=%.6f moves=%d\n",
               sweep, cur_hpwl, overflow_ratio, overflow_threshold, moves_this_sweep);

        if (cur_hpwl < best_hpwl) best_hpwl = cur_hpwl;

        // Check convergence: relax if stuck
        float_t hpwl_ratio = (sweep > 0) ? cur_hpwl / best_hpwl : 1.0f;
        if (false && (moves_this_sweep < 10 || hpwl_ratio > cwtr_cfg.hpwl_improve_min)) {
            if (overflow_threshold < 1.0f) {
                overflow_threshold = std::min(1.0f,
                    overflow_threshold * cwtr_cfg.overflow_relax);
                printf("  → Relaxing overflow threshold to %.4f\n", overflow_threshold);
            } else {
                printf("  → Converged (threshold at max, no improvement)\n");
                break;
            }
        } else if (max_overflow < overflow_threshold * 0.8f) {
            // Tighten constraint if we're doing better than threshold
            overflow_threshold = std::max(0.01f, max_overflow * 1.1f);
        }

        // Swap operations when stuck
        if (cwtr_cfg.swap_attempts > 0 && moves_this_sweep < 5 && sweep > 3) {
            printf("  → Attempting swap operations...\n");
            int swaps_done = 0;
            std::uniform_int_distribution<int> dist(0, Nm-1);

            for (int s = 0; s < cwtr_cfg.swap_attempts && swaps_done < 50; ++s) {
                int a = movable_ids[dist(rng)];
                int b = movable_ids[dist(rng)];
                if (a == b) continue;

                float_t ax = cells[a].x, ay = cells[a].y;
                float_t bx = cells[b].x, by = cells[b].y;

                // Tentative swap: a→b_pos, b→a_pos
                float_t ov_a_at_b = compute_overflow_after_move(
                    a, bx, by, cells, grid, bin_area_sum, current_excess_area) /
                    grid.available_area;
                float_t ov_b_at_a = compute_overflow_after_move(
                    b, ax, ay, cells, grid, bin_area_sum, current_excess_area) /
                    grid.available_area;

                if (ov_a_at_b <= overflow_threshold && ov_b_at_a <= overflow_threshold) {
                    float_t hpwl_delta = compute_hpwl_delta(a, bx, by, cells, nets, adj)
                                       + compute_hpwl_delta(b, ax, ay, cells, nets, adj);

                    if (hpwl_delta < -1.0f) {
                        // Execute swap
                        float_t hw_a = cells[a].width * 0.5f, hh_a = cells[a].height * 0.5f;
                        float_t hw_b = cells[b].width * 0.5f, hh_b = cells[b].height * 0.5f;

                        // Remove both
                        auto update_bins_cell = [&](int /*cid*/, float_t cx, float_t cy,
                                                      float_t hw, float_t hh, float_t sign) {
                            float_t cl = cx - hw, cr = cx + hw;
                            float_t cb = cy - hh, ct = cy + hh;
                            int bx_s = std::max(0,(int)((cl-grid.chip_xl)/grid.bin_w));
                            int bx_e = std::min(grid.nx-1,(int)((cr-grid.chip_xl)/grid.bin_w));
                            int by_s = std::max(0,(int)((cb-grid.chip_yl)/grid.bin_h));
                            int by_e = std::min(grid.ny-1,(int)((ct-grid.chip_yl)/grid.bin_h));
                            for (int by = by_s; by <= by_e; ++by) {
                                float_t byl = grid.chip_yl+by*grid.bin_h;
                                float_t oy = std::max(0.0f,std::min(ct,byl+grid.bin_h)-std::max(cb,byl));
                                for (int bx = bx_s; bx <= bx_e; ++bx) {
                                    float_t bxl = grid.chip_xl+bx*grid.bin_w;
                                    float_t ox = std::max(0.0f,std::min(cr,bxl+grid.bin_w)-std::max(cl,bxl));
                                    bin_area_sum[by*grid.nx+bx] += sign*ox*oy;
                                }
                            }
                        };
                        update_bins_cell(a, ax, ay, hw_a, hh_a, -1.0f);
                        update_bins_cell(b, bx, by, hw_b, hh_b, -1.0f);
                        cells[a].x = bx; cells[a].y = by;
                        cells[b].x = ax; cells[b].y = ay;
                        update_bins_cell(a, bx, by, hw_a, hh_a, +1.0f);
                        update_bins_cell(b, ax, ay, hw_b, hh_b, +1.0f);

                        ++swaps_done;
                        ++total_moves;
                    }
                }
            }
            printf("  → Swaps accepted: %d\n", swaps_done);
        }
    }

    csv_log.close();
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();

    float_t final_hpwl = compute_hpwl_only(cells, nets);
    float_t final_overflow = density_evaluate(grid, cells);
    if (final_overflow > cwtr_cfg.overflow_limit + 1.0e-7f) {
        printf("[CWTR] Rolling back: recomputed overflow %.6f exceeds %.6f\n",
               final_overflow, cwtr_cfg.overflow_limit);
        cells = initial_cells;
        final_hpwl = initial_hpwl;
        final_overflow = density_evaluate(grid, cells);
        total_moves = 0;
    }

    printf("\n[CWTR] Complete: %d moves, %.2f sec\n", total_moves, elapsed);
    printf("[CWTR] HPWL: %.1f → %.1f  (%.1f%% reduction)\n",
           initial_hpwl, final_hpwl, (initial_hpwl-final_hpwl)/initial_hpwl*100);
    printf("[CWTR] Overflow: %.6f\n", final_overflow);

    return { final_hpwl, final_overflow, total_moves, elapsed };
}
