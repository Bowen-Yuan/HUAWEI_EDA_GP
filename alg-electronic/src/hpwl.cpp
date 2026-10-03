#include "types.h"
#include <cmath>
#include <cstdio>
#include <omp.h>

namespace {
int hpwl_gradient_net_degree_limit = 0;
}

void set_hpwl_gradient_net_degree_limit(int degree_limit) {
    hpwl_gradient_net_degree_limit = std::max(0, degree_limit);
}

// Compute exact HPWL and its subgradient (no smoothing whatsoever).
// For each net: HPWL = max_x - min_x + max_y - min_y
// Subgradient:
//   ∂/∂x_i = +1 (sole max), -1 (sole min), or proportional share for ties
//   ∂/∂y_i = +1 (sole max_y), -1 (sole min_y), or proportional share for ties
//
// Performance optimizations:
//   - Flat arrays for positions (SoA layout for SIMD-friendly access)
//   - OpenMP parallel over nets with per-thread gradient accumulation
//   - Branch-minimized inner loops

float_t compute_hpwl_and_gradient(
    const std::vector<Cell>& cells,
    const std::vector<Net>& nets,
    std::vector<float_t>& grad_x,
    std::vector<float_t>& grad_y)
{
    const int N = (int)cells.size();
    const int M = (int)nets.size();
    static bool reported_parallel_shape = false;
    if (!reported_parallel_shape) {
        std::printf("[ExactHPWL] nodes=%d nets=%d omp_max_threads=%d\n",
                    N, M, omp_get_max_threads());
        reported_parallel_shape = true;
    }

    // Reset gradients (only for movable cells)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        if (!cells[i].is_terminal) {
            grad_x[i] = 0.0f;
            grad_y[i] = 0.0f;
        }
    }

    float_t total_hpwl = 0.0f;

    #pragma omp parallel
    {
        // Reuse each worker's buffers across iterations. Reallocating two
        // N-sized arrays per OpenMP worker on every exact-HPWL evaluation can
        // fragment the Windows heap during long placement runs.
        static thread_local std::vector<float_t> local_gx;
        static thread_local std::vector<float_t> local_gy;
        local_gx.assign(N, 0.0f);
        local_gy.assign(N, 0.0f);
        float_t local_hpwl = 0.0f;

        #pragma omp for schedule(dynamic, 256)
        for (int e = 0; e < M; ++e) {
            const auto& net = nets[e];
            const int npins = (int)net.cell_ids.size();
            if (npins < 2) continue;

            // --- Find min/max x and y with pin indices ---
            int first_id = net.cell_ids[0];
            float_t min_x = cells[first_id].x;
            float_t max_x = min_x;
            float_t min_y = cells[first_id].y;
            float_t max_y = min_y;

            int n_min_x = 1, n_max_x = 1;
            int n_min_y = 1, n_max_y = 1;

            // Single pass to find extrema
            for (int p = 1; p < npins; ++p) {
                int cid = net.cell_ids[p];
                float_t cx = cells[cid].x;
                float_t cy = cells[cid].y;

                if (cx < min_x)      { min_x = cx; n_min_x = 1; }
                else if (cx == min_x) { ++n_min_x; }

                if (cx > max_x)      { max_x = cx; n_max_x = 1; }
                else if (cx == max_x) { ++n_max_x; }

                if (cy < min_y)      { min_y = cy; n_min_y = 1; }
                else if (cy == min_y) { ++n_min_y; }

                if (cy > max_y)      { max_y = cy; n_max_y = 1; }
                else if (cy == max_y) { ++n_max_y; }
            }

            float_t hpwl_e = (max_x - min_x) + (max_y - min_y);
            local_hpwl += hpwl_e;

            // DREAMPlace excludes very high-degree nets from the optimization
            // direction. They remain in local_hpwl, so reporting and feasible
            // state selection still use the exact all-net HPWL objective.
            if (hpwl_gradient_net_degree_limit > 0 &&
                npins >= hpwl_gradient_net_degree_limit)
                continue;

            // --- Accumulate subgradients ---
            // Gradient shares: +1/n_max for max pins, -1/n_min for min pins
            const float_t gx_max_val = 1.0f / n_max_x;
            const float_t gx_min_val = -1.0f / n_min_x;
            const float_t gy_max_val = 1.0f / n_max_y;
            const float_t gy_min_val = -1.0f / n_min_y;

            // Second pass: assign gradients
            for (int p = 0; p < npins; ++p) {
                int cid = net.cell_ids[p];
                if (cells[cid].is_terminal) continue;

                float_t cx = cells[cid].x;
                float_t cy = cells[cid].y;

                // X gradient
                float_t gx = 0.0f;
                if (cx == max_x) gx += gx_max_val;
                if (cx == min_x) gx += gx_min_val;
                local_gx[cid] += gx;

                // Y gradient
                float_t gy = 0.0f;
                if (cy == max_y) gy += gy_max_val;
                if (cy == min_y) gy += gy_min_val;
                local_gy[cid] += gy;
            }
        }

        // Merge thread-local gradients into global arrays
        #pragma omp critical
        {
            total_hpwl += local_hpwl;
            for (int i = 0; i < N; ++i) {
                if (!cells[i].is_terminal) {
                    grad_x[i] += local_gx[i];
                    grad_y[i] += local_gy[i];
                }
            }
        }
    }

    return total_hpwl;
}

// Fast HPWL evaluation only (no gradient) — used for reporting
float_t compute_hpwl_only(const std::vector<Cell>& cells, const std::vector<Net>& nets)
{
    const int M = (int)nets.size();
    float_t total = 0.0f;

    #pragma omp parallel for schedule(dynamic, 256) reduction(+:total)
    for (int e = 0; e < M; ++e) {
        const auto& net = nets[e];
        if (net.cell_ids.size() < 2) continue;

        float_t min_x = cells[net.cell_ids[0]].x;
        float_t max_x = min_x;
        float_t min_y = cells[net.cell_ids[0]].y;
        float_t max_y = min_y;

        for (size_t p = 1; p < net.cell_ids.size(); ++p) {
            float_t cx = cells[net.cell_ids[p]].x;
            float_t cy = cells[net.cell_ids[p]].y;
            if (cx < min_x) min_x = cx;
            if (cx > max_x) max_x = cx;
            if (cy < min_y) min_y = cy;
            if (cy > max_y) max_y = cy;
        }
        total += (max_x - min_x) + (max_y - min_y);
    }
    return total;
}
