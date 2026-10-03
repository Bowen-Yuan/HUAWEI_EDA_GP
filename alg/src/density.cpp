#include "types.h"
#include <cmath>
#include <omp.h>
#include <algorithm>

// ============================================================================
// Non-smooth density penalty: D(x,y) = sum_b max(0, rho_b - rho_target)^2
// With exact overlap computation -- no Gaussian or any other smoothing.
// ============================================================================

// Runtime-configurable density gradient scale (default: 0.02 for balanced WL-density)
// Set to higher values (0.5-2.0) for aggressive spreading phases
float_t g_density_grad_scale = 0.02f;
bool   g_short_edge_only = false;   // Only keep the dominant gradient direction per cell

void set_density_grad_scale(float_t s) { g_density_grad_scale = s; }
float_t get_density_grad_scale() { return g_density_grad_scale; }
void set_short_edge_only(bool v) { g_short_edge_only = v; }

// Initialize bin grid based on chip dimensions and cell sizes
void density_init(BinGrid& grid,
                  float_t chip_xl, float_t chip_yl,
                  float_t chip_xh, float_t chip_yh,
                  const std::vector<Cell>& cells,
                  float_t target_density)
{
    grid.chip_xl = chip_xl;
    grid.chip_yl = chip_yl;
    grid.chip_xh = chip_xh;
    grid.chip_yh = chip_yh;
    grid.target_density = target_density;

    float_t chip_w = chip_xh - chip_xl;
    float_t chip_h = chip_yh - chip_yl;

    // Compute average cell dimensions for bin sizing
    float_t avg_w = 0, avg_h = 0;
    int n_movable = 0;
    for (const auto& c : cells) {
        if (!c.is_terminal) {
            avg_w += c.width;
            avg_h += c.height;
            ++n_movable;
        }
    }
    if (n_movable == 0) { n_movable = 1; avg_w = chip_w / 100; avg_h = chip_h / 100; }
    else { avg_w /= n_movable; avg_h /= n_movable; }

    // Bin size: roughly 10-20x average cell size for good density resolution
    float_t bin_w_target = std::max(avg_w * 15.0f, chip_w / 200.0f);
    float_t bin_h_target = std::max(avg_h * 15.0f, chip_h / 200.0f);

    grid.nx = std::max(8, (int)(chip_w / bin_w_target));
    grid.ny = std::max(8, (int)(chip_h / bin_h_target));
    grid.bin_w = chip_w / grid.nx;
    grid.bin_h = chip_h / grid.ny;

    int total_bins = grid.nx * grid.ny;
    grid.density.assign(total_bins, 0.0f);
    grid.overflow.assign(total_bins, 0.0f);
    grid.total_overflow = 0.0f;

    printf("[Density] Grid: %d x %d = %d bins, bin size %.1f x %.1f, target_density=%.4f\n",
           grid.nx, grid.ny, total_bins, grid.bin_w, grid.bin_h, target_density);
}

// Compute density, overflow penalty, and subgradient (exact, no smoothing)
// Returns total density penalty D = Σ max(0, ρ - ρ_tgt)²
float_t density_compute_and_gradient(
    const BinGrid& grid,
    const std::vector<Cell>& cells,
    std::vector<float_t>& grad_x,
    std::vector<float_t>& grad_y,
    std::vector<float_t>& bin_density)
{
    const int N = (int)cells.size();
    const int nbins = grid.nx * grid.ny;
    const float_t bin_area = grid.bin_w * grid.bin_h;
    const float_t inv_bin_area = 1.0f / bin_area;
    const float_t rho_target = grid.target_density;

    // Zero out bin densities
    bin_density.assign(nbins, 0.0f);

    // ---- Phase 1: Compute bin densities ----
    // Each thread accumulates into a private bin array to avoid atomics
    #pragma omp parallel
    {
        std::vector<float_t> local_bins(nbins, 0.0f);

        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < N; ++i) {
            if (cells[i].is_terminal) continue;

            float_t cx = cells[i].x;
            float_t cy = cells[i].y;
            float_t hw = cells[i].width * 0.5f;
            float_t hh = cells[i].height * 0.5f;

            // Cell bounding box
            float_t cl = cx - hw, cr = cx + hw;
            float_t cb = cy - hh, ct = cy + hh;

            // Bin range that this cell overlaps
            int bx_start = (int)((cl - grid.chip_xl) / grid.bin_w);
            int bx_end   = (int)((cr - grid.chip_xl) / grid.bin_w);
            int by_start = (int)((cb - grid.chip_yl) / grid.bin_h);
            int by_end   = (int)((ct - grid.chip_yl) / grid.bin_h);

            bx_start = std::max(0, bx_start);
            bx_end   = std::min(grid.nx - 1, bx_end);
            by_start = std::max(0, by_start);
            by_end   = std::min(grid.ny - 1, by_end);

            for (int by = by_start; by <= by_end; ++by) {
                float_t bin_yl = grid.chip_yl + by * grid.bin_h;
                float_t bin_yh = bin_yl + grid.bin_h;
                float_t oy = std::max(0.0f, std::min(ct, bin_yh) - std::max(cb, bin_yl));

                for (int bx = bx_start; bx <= bx_end; ++bx) {
                    float_t bin_xl = grid.chip_xl + bx * grid.bin_w;
                    float_t bin_xh = bin_xl + grid.bin_w;
                    float_t ox = std::max(0.0f, std::min(cr, bin_xh) - std::max(cl, bin_xl));

                    int bidx = by * grid.nx + bx;
                    local_bins[bidx] += ox * oy;
                }
            }
        }

        #pragma omp critical
        {
            for (int b = 0; b < nbins; ++b)
                bin_density[b] += local_bins[b];
        }
    }

    // ---- Phase 2: Compute density penalty and per-bin overflow ----
    // D = Σ max(0, ρ-ρ_target)²  [dimensionless, ρ is density ratio]
    // ∂D/∂x_i = 2*Σ_b max(0,ρ-ρ_target) * (1/bin_area) * (±1) * oy
    // This is O(2*over*oy/bin_area) ≈ O(1e-5) — far too small vs HPWL gradient O(1).
    // Scale by bin_area to bring to comparable magnitude => og = 2*over (unitless).
    // A further global scale factor (default ~0.2) allows tuning density vs WL balance.
    float_t total_penalty = 0.0f;
    std::vector<float_t> overflow_grad(nbins, 0.0f);
    float_t grad_scale = g_density_grad_scale;  // local copy for thread safety

    #pragma omp parallel for schedule(static) reduction(+:total_penalty)
    for (int b = 0; b < nbins; ++b) {
        float_t rho = bin_density[b] * inv_bin_area;
        float_t over = rho - rho_target;
        if (over > 0.0f) {
            total_penalty += over * over;
            overflow_grad[b] = 2.0f * over * grad_scale;
        } else {
            overflow_grad[b] = 0.0f;
        }
    }

    // ---- Phase 3: Back-propagate bin gradients to cell positions ----
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < N; ++i) {
        if (cells[i].is_terminal) continue;

        float_t cx = cells[i].x;
        float_t cy = cells[i].y;
        float_t hw = cells[i].width * 0.5f;
        float_t hh = cells[i].height * 0.5f;

        float_t cl = cx - hw, cr = cx + hw;
        float_t cb = cy - hh, ct = cy + hh;

        int bx_start = (int)((cl - grid.chip_xl) / grid.bin_w);
        int bx_end   = (int)((cr - grid.chip_xl) / grid.bin_w);
        int by_start = (int)((cb - grid.chip_yl) / grid.bin_h);
        int by_end   = (int)((ct - grid.chip_yl) / grid.bin_h);

        bx_start = std::max(0, bx_start);
        bx_end   = std::min(grid.nx - 1, bx_end);
        by_start = std::max(0, by_start);
        by_end   = std::min(grid.ny - 1, by_end);

        float_t dgx = 0.0f, dgy = 0.0f;

        for (int by = by_start; by <= by_end; ++by) {
            float_t bin_yl = grid.chip_yl + by * grid.bin_h;
            float_t bin_yh = bin_yl + grid.bin_h;
            float_t oy = std::max(0.0f, std::min(ct, bin_yh) - std::max(cb, bin_yl));

            // Derivative of oy w.r.t. cy
            // oy = min(ct, by_t) - max(cb, by_b)
            float_t doy_dcy = 0.0f;
            if (oy > 0.0f) {
                int cbin = (cb >= bin_yl && cb < bin_yh) ? 1 : 0;  // bottom edge in
                int ctin = (ct > bin_yl && ct <= bin_yh) ? 1 : 0;  // top edge in
                if (cbin && !ctin)  doy_dcy = -1.0f;  // bottom inside, top outside
                if (!cbin && ctin) doy_dcy = +1.0f;  // top inside, bottom outside
            }

            for (int bx = bx_start; bx <= bx_end; ++bx) {
                int bidx = by * grid.nx + bx;
                float_t og = overflow_grad[bidx];
                if (og == 0.0f) continue;

                float_t bin_xl = grid.chip_xl + bx * grid.bin_w;
                float_t bin_xh = bin_xl + grid.bin_w;
                float_t ox = std::max(0.0f, std::min(cr, bin_xh) - std::max(cl, bin_xl));

                if (ox <= 0.0f || oy <= 0.0f) continue;

                // Derivative of ox w.r.t. cx
                float_t dox_dcx = 0.0f;
                int clin = (cl >= bin_xl && cl < bin_xh) ? 1 : 0;
                int crin = (cr > bin_xl && cr <= bin_xh) ? 1 : 0;
                if (clin && !crin)  dox_dcx = -1.0f;  // left inside, right outside
                if (!clin && crin) dox_dcx = +1.0f;  // right inside, left outside

                // Chain rule: ∂D/∂cx = og * ∂(area)/∂cx
                // area = ox * oy
                // ∂area/∂cx = dox_dcx * oy
                dgx += og * dox_dcx * oy;

                // ∂area/∂cy = ox * doy_dcy
                dgy += og * ox * doy_dcy;
            }
        }

        // Short-edge-only: zero out less effective direction
        // dgx ∝ overlap_y, dgy ∝ overlap_x.
        // If |dgx| > |dgy|, x is the short-edge direction → keep only x.
        if (g_short_edge_only) {
            if (std::fabs(dgx) > std::fabs(dgy))
                dgy = 0.0f;
            else
                dgx = 0.0f;
        }

        grad_x[i] += dgx;
        grad_y[i] += dgy;
    }

    return total_penalty;
}

// Fast density evaluation (no gradient) for reporting
float_t density_evaluate(const BinGrid& grid, const std::vector<Cell>& cells)
{
    const int nbins = grid.nx * grid.ny;
    const float_t bin_area = grid.bin_w * grid.bin_h;
    const float_t inv_bin_area = 1.0f / bin_area;

    std::vector<float_t> bins(nbins, 0.0f);

    #pragma omp parallel
    {
        std::vector<float_t> local(nbins, 0.0f);
        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < (int)cells.size(); ++i) {
            if (cells[i].is_terminal) continue;

            float_t cx = cells[i].x, cy = cells[i].y;
            float_t hw = cells[i].width * 0.5f, hh = cells[i].height * 0.5f;
            float_t cl = cx - hw, cr = cx + hw;
            float_t cb = cy - hh, ct = cy + hh;

            int bx_s = std::max(0, (int)((cl - grid.chip_xl) / grid.bin_w));
            int bx_e = std::min(grid.nx-1, (int)((cr - grid.chip_xl) / grid.bin_w));
            int by_s = std::max(0, (int)((cb - grid.chip_yl) / grid.bin_h));
            int by_e = std::min(grid.ny-1, (int)((ct - grid.chip_yl) / grid.bin_h));

            for (int by = by_s; by <= by_e; ++by) {
                float_t oy = std::max(0.0f,
                    std::min(ct, grid.chip_yl+(by+1)*grid.bin_h) -
                    std::max(cb, grid.chip_yl+by*grid.bin_h));
                for (int bx = bx_s; bx <= bx_e; ++bx) {
                    float_t ox = std::max(0.0f,
                        std::min(cr, grid.chip_xl+(bx+1)*grid.bin_w) -
                        std::max(cl, grid.chip_xl+bx*grid.bin_w));
                    local[by*grid.nx + bx] += ox * oy;
                }
            }
        }
        #pragma omp critical
        { for (int b=0; b<nbins; ++b) bins[b] += local[b]; }
    }

    float_t total_over = 0.0f;
    float_t max_over = 0.0f;
    for (int b = 0; b < nbins; ++b) {
        float_t rho = bins[b] * inv_bin_area;
        float_t over = rho - grid.target_density;
        if (over > 0.0f) {
            total_over += over;
            if (over > max_over) max_over = over;
        }
    }
    printf("[Density] avg_overflow=%.6f, max_overflow=%.6f\n",
           total_over / nbins, max_over);
    return total_over / nbins;
}
