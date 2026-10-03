#include "types.h"
#include "bayesian_lambda.h"
#include "baseline_lambda.h"
#include "trajectory_lambda.h"
#include "visualization.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <fstream>
#include <limits>
#include <utility>

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
    double fixed_blockage_area = 0.0;
    for (const Cell& cell : cells) {
        if (!cell.is_terminal) continue;
        const double left = std::max((double)chip_xl,
                                     (double)cell.x - 0.5 * cell.width);
        const double right = std::min((double)chip_xh,
                                      (double)cell.x + 0.5 * cell.width);
        const double bottom = std::max((double)chip_yl,
                                       (double)cell.y - 0.5 * cell.height);
        const double top = std::min((double)chip_yh,
                                    (double)cell.y + 0.5 * cell.height);
        if (right > left && top > bottom)
            fixed_blockage_area += (right - left) * (top - bottom);
    }
    fixed_blockage_area = std::min(fixed_blockage_area, chip_area * 0.99);
    const double available_area = chip_area - fixed_blockage_area;

    // Match the project's baseline metric. Fixed cells are represented by
    // target-density blockage charge, which is algebraically equivalent to
    // subtracting their area from each bin's movable capacity.
    float_t target_density = cfg.target_density_override > 0.0f
        ? cfg.target_density_override
        : 0.80f;

    printf("[Solver] movable_density=%.4f fixed_blockage=%.4f target=%.4f%s\n",
           total_cell_area / chip_area, fixed_blockage_area / chip_area,
           target_density, cfg.target_density_override > 0.0f ? " (override)" : "");

    std::vector<float_t> grad_x(N, 0.0f);
    std::vector<float_t> grad_y(N, 0.0f);
    std::vector<float_t> density_grad_x(N, 0.0f);
    std::vector<float_t> density_grad_y(N, 0.0f);
    std::vector<float_t> mom_x(N, 0.0f);
    std::vector<float_t> mom_y(N, 0.0f);
    std::vector<float_t> second_x(N, 0.0f);
    std::vector<float_t> second_y(N, 0.0f);

    struct HpwlBundleCut {
        float_t hpwl = 0.0f;
        std::vector<float_t> x, y, gx, gy;
    };
    struct HpwlBundleStats {
        bool updated = false;
        int cuts = 0;
        float_t tau = 0.0f;
        float_t newest_weight = 1.0f;
        float_t model_error = 0.0f;
        float_t norm_ratio = 1.0f;
        float_t entropy = 0.0f;
        float_t trigger_cosine = 1.0f;
        int cut_age = 0;
        float_t effective_current_mix = 1.0f;
    };
    std::vector<HpwlBundleCut> hpwl_bundle;
    std::vector<float_t> bundle_grad_x(Nm, 0.0f);
    std::vector<float_t> bundle_grad_y(Nm, 0.0f);
    bool bundle_direction_ready = false;
    int bundle_last_cut_iter = -1000000;
    std::ofstream bundle_log;
    if (cfg.use_hpwl_bundle) {
        bundle_log.open("bundle_strategy.csv");
        bundle_log << "global_iteration,fine_iteration,cuts,tau,newest_weight,"
                      "model_error,norm_ratio,entropy,current_mix,trigger_cosine,cut_age\n";
    }

    // ---- Bin grids ----
    BinGrid grid_fine, grid_medium, grid_coarse;
    density_init(grid_fine, chip_xl, chip_yl, chip_xh, chip_yh, cells, target_density);
    grid_fine.available_area = (float_t)available_area;
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

    // DREAMPlace initializes the density weight from the gradient-norm ratio.
    // Keep the user-facing lambda schedule dimensionless and map it to the
    // physical electric-energy coefficient once, at the first density step.
    float_t density_weight_scale = -1.0f;
    auto combine_density_gradient = [&](float_t lambda_control) -> float_t {
        const bool first_scale_for_grid = density_weight_scale < 0.0f;
        if (density_weight_scale < 0.0f) {
            double wl_norm = 0.0;
            double density_norm = 0.0;
            for (int idx = 0; idx < Nm; ++idx) {
                const int i = movable_ids[idx];
                wl_norm += std::fabs(grad_x[i]) + std::fabs(grad_y[i]);
                density_norm += std::fabs(density_grad_x[i]) +
                                std::fabs(density_grad_y[i]);
            }
            density_weight_scale = density_norm > 1e-20
                ? (float_t)(wl_norm / density_norm) : 0.0f;
            if (first_scale_for_grid)
                printf("[ElectricDensity] gradient normalization scale=%.6e\n",
                       density_weight_scale);
        }
        const float_t effective_lambda = lambda_control * density_weight_scale;
        for (int idx = 0; idx < Nm; ++idx) {
            const int i = movable_ids[idx];
            grad_x[i] += effective_lambda * density_grad_x[i];
            grad_y[i] += effective_lambda * density_grad_y[i];
        }
        return effective_lambda;
    };

    auto project_hpwl_from_density_ascent = [&]() -> float_t {
        double dot = 0.0;
        double density_norm2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            const int i = movable_ids[idx];
            dot += (double)grad_x[i] * density_grad_x[i] +
                   (double)grad_y[i] * density_grad_y[i];
            density_norm2 += (double)density_grad_x[i] * density_grad_x[i] +
                             (double)density_grad_y[i] * density_grad_y[i];
        }
        if (dot >= 0.0 || density_norm2 <= 1.0e-30)
            return 0.0f;

        const float_t coefficient = (float_t)(dot / density_norm2);
        for (int idx = 0; idx < Nm; ++idx) {
            const int i = movable_ids[idx];
            grad_x[i] -= coefficient * density_grad_x[i];
            grad_y[i] -= coefficient * density_grad_y[i];
        }
        return -coefficient;
    };

    auto apply_hpwl_bundle = [&](int fine_iter, float_t current_hpwl,
                                 float_t lr, float_t hpwl_rms) {
        HpwlBundleStats stats;
        if (!cfg.use_hpwl_bundle) return stats;

        float_t trigger_cosine = 1.0f;
        const int cut_age = fine_iter - bundle_last_cut_iter;
        if (!hpwl_bundle.empty()) {
            const HpwlBundleCut& newest = hpwl_bundle.back();
            double dot = 0.0;
            double current_norm2 = 0.0;
            double newest_norm2 = 0.0;
            for (int idx = 0; idx < Nm; ++idx) {
                const int i = movable_ids[idx];
                dot += (double)grad_x[i] * newest.gx[idx] +
                       (double)grad_y[i] * newest.gy[idx];
                current_norm2 += (double)grad_x[i] * grad_x[i] +
                                 (double)grad_y[i] * grad_y[i];
                newest_norm2 += (double)newest.gx[idx] * newest.gx[idx] +
                                (double)newest.gy[idx] * newest.gy[idx];
            }
            const double denominator = std::sqrt(current_norm2 * newest_norm2);
            trigger_cosine = denominator > 1.0e-20
                ? (float_t)(dot / denominator) : 1.0f;
        }
        bool insert_cut = fine_iter % cfg.hpwl_bundle_interval == 0;
        if (cfg.hpwl_bundle_adaptive_cuts) {
            insert_cut = hpwl_bundle.empty() || cut_age >= cfg.hpwl_bundle_interval;
            if (!hpwl_bundle.empty() && cut_age >= cfg.hpwl_bundle_min_interval) {
                if (trigger_cosine <= cfg.hpwl_bundle_cosine_trigger)
                    insert_cut = true;
            }
        }
        if (insert_cut) {
            HpwlBundleCut cut;
            cut.hpwl = current_hpwl;
            cut.x.resize(Nm);
            cut.y.resize(Nm);
            cut.gx.resize(Nm);
            cut.gy.resize(Nm);
            for (int idx = 0; idx < Nm; ++idx) {
                const int i = movable_ids[idx];
                cut.x[idx] = cells[i].x;
                cut.y[idx] = cells[i].y;
                cut.gx[idx] = grad_x[i];
                cut.gy[idx] = grad_y[i];
            }
            if ((int)hpwl_bundle.size() >= cfg.hpwl_bundle_size)
                hpwl_bundle.erase(hpwl_bundle.begin());
            hpwl_bundle.push_back(std::move(cut));
            bundle_last_cut_iter = fine_iter;

            const int K = (int)hpwl_bundle.size();
            std::vector<double> beta(K, 0.0);
            std::vector<double> gram(K * K, 0.0);
            for (int j = 0; j < K; ++j) {
                const HpwlBundleCut& cj = hpwl_bundle[j];
                double affine_shift = 0.0;
                for (int idx = 0; idx < Nm; ++idx) {
                    const int i = movable_ids[idx];
                    affine_shift += (double)cj.gx[idx] * (cells[i].x - cj.x[idx]) +
                                    (double)cj.gy[idx] * (cells[i].y - cj.y[idx]);
                }
                beta[j] = (double)cj.hpwl + affine_shift;
                for (int k = 0; k <= j; ++k) {
                    const HpwlBundleCut& ck = hpwl_bundle[k];
                    double dot = 0.0;
                    for (int idx = 0; idx < Nm; ++idx)
                        dot += (double)cj.gx[idx] * ck.gx[idx] +
                               (double)cj.gy[idx] * ck.gy[idx];
                    gram[j * K + k] = dot;
                    gram[k * K + j] = dot;
                }
            }

            const double tau = cfg.hpwl_bundle_prox_scale * lr /
                               std::max((double)hpwl_rms, 1.0e-6);
            std::vector<double> alpha(K, 0.0);
            alpha[K - 1] = 1.0;
            for (int fw_iter = 0; fw_iter < 20; ++fw_iter) {
                std::vector<double> dual_gradient(K, 0.0);
                int best = 0;
                for (int j = 0; j < K; ++j) {
                    double q_alpha = 0.0;
                    for (int k = 0; k < K; ++k)
                        q_alpha += gram[j * K + k] * alpha[k];
                    dual_gradient[j] = beta[j] - tau * q_alpha;
                    if (dual_gradient[j] > dual_gradient[best]) best = j;
                }
                std::vector<double> direction(K, 0.0);
                for (int j = 0; j < K; ++j) direction[j] = -alpha[j];
                direction[best] += 1.0;
                double numerator = 0.0;
                double curvature = 0.0;
                for (int j = 0; j < K; ++j) {
                    numerator += dual_gradient[j] * direction[j];
                    for (int k = 0; k < K; ++k)
                        curvature += direction[j] * gram[j * K + k] * direction[k];
                }
                if (numerator <= 1.0e-7) break;
                const double gamma = curvature > 1.0e-20
                    ? std::min(1.0, numerator / (tau * curvature)) : 1.0;
                for (int j = 0; j < K; ++j)
                    alpha[j] += gamma * direction[j];
            }

            std::fill(bundle_grad_x.begin(), bundle_grad_x.end(), 0.0f);
            std::fill(bundle_grad_y.begin(), bundle_grad_y.end(), 0.0f);
            double aggregate_beta = 0.0;
            double aggregate_norm2 = 0.0;
            double current_norm2 = 0.0;
            double entropy = 0.0;
            for (int j = 0; j < K; ++j) {
                aggregate_beta += alpha[j] * beta[j];
                if (alpha[j] > 1.0e-15) entropy -= alpha[j] * std::log(alpha[j]);
                for (int idx = 0; idx < Nm; ++idx) {
                    bundle_grad_x[idx] += (float_t)alpha[j] * hpwl_bundle[j].gx[idx];
                    bundle_grad_y[idx] += (float_t)alpha[j] * hpwl_bundle[j].gy[idx];
                }
            }
            for (int idx = 0; idx < Nm; ++idx) {
                const int i = movable_ids[idx];
                aggregate_norm2 += (double)bundle_grad_x[idx] * bundle_grad_x[idx] +
                                   (double)bundle_grad_y[idx] * bundle_grad_y[idx];
                current_norm2 += (double)grad_x[i] * grad_x[i] +
                                 (double)grad_y[i] * grad_y[i];
            }
            bundle_direction_ready = true;
            stats.updated = true;
            stats.cuts = K;
            stats.tau = (float_t)tau;
            stats.newest_weight = (float_t)alpha[K - 1];
            stats.model_error = (float_t)std::max(0.0,
                (double)current_hpwl - aggregate_beta);
            stats.norm_ratio = current_norm2 > 1.0e-20
                ? (float_t)std::sqrt(aggregate_norm2 / current_norm2) : 1.0f;
            stats.entropy = (float_t)entropy;
            stats.trigger_cosine = trigger_cosine;
            stats.cut_age = cut_age;
        }

        float_t effective_current_mix = cfg.hpwl_bundle_current_mix;
        if (cfg.hpwl_bundle_extrema_mix) {
            float_t gain = cfg.hpwl_bundle_mix_gain;
            if (cfg.hpwl_bundle_mix_decay_steps > 0) {
                gain *= std::max(0.0f, 1.0f - (float_t)fine_iter /
                    cfg.hpwl_bundle_mix_decay_steps);
            }
            effective_current_mix = std::max(0.20f, std::min(0.80f,
                cfg.hpwl_bundle_current_mix +
                gain * trigger_cosine));
        }
        if (bundle_direction_ready) {
            const float_t current_mix = effective_current_mix;
            const float_t bundle_mix = 1.0f - current_mix;
            for (int idx = 0; idx < Nm; ++idx) {
                const int i = movable_ids[idx];
                grad_x[i] = current_mix * grad_x[i] + bundle_mix * bundle_grad_x[idx];
                grad_y[i] = current_mix * grad_y[i] + bundle_mix * bundle_grad_y[idx];
            }
        }
        stats.effective_current_mix = effective_current_mix;
        return stats;
    };

    float_t chip_w = chip_xh - chip_xl, chip_h = chip_yh - chip_yl;
    float_t hpwl = 0.0f, dpen = 0.0f;
    float_t hpwl_best = 1e30f;
    int total_iters = 0;
    const float_t RMS_EPS = 1e-6f;
    int adam_age = 0;

    auto reset_density_optimizer = [&]() {
        if (!cfg.use_adam) return;
        std::fill(mom_x.begin(), mom_x.end(), 0.0f);
        std::fill(mom_y.begin(), mom_y.end(), 0.0f);
        std::fill(second_x.begin(), second_x.end(), 0.0f);
        std::fill(second_y.begin(), second_y.end(), 0.0f);
        adam_age = 0;
    };

    auto apply_density_step = [&](float_t g_rms, float_t lr,
                                  float_t momentum, bool allow_margin) {
        if (cfg.use_adam) ++adam_age;
        const float_t beta1 = 0.80f;
        const float_t beta2 = 0.99f;
        const float_t b1_correction = cfg.use_adam
            ? std::max(1.0e-8f, 1.0f - std::pow(beta1, (float_t)adam_age)) : 1.0f;
        const float_t b2_correction = cfg.use_adam
            ? std::max(1.0e-8f, 1.0f - std::pow(beta2, (float_t)adam_age)) : 1.0f;

        for (int idx = 0; idx < Nm; ++idx) {
            const int i = movable_ids[idx];
            float_t gx = grad_x[i] / g_rms;
            float_t gy = grad_y[i] / g_rms;
            gx = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gx));
            gy = std::max(-cfg.grad_clip, std::min(cfg.grad_clip, gy));

            float_t step_x = 0.0f;
            float_t step_y = 0.0f;
            if (cfg.use_adam) {
                mom_x[i] = beta1 * mom_x[i] + (1.0f - beta1) * gx;
                mom_y[i] = beta1 * mom_y[i] + (1.0f - beta1) * gy;
                second_x[i] = beta2 * second_x[i] + (1.0f - beta2) * gx * gx;
                second_y[i] = beta2 * second_y[i] + (1.0f - beta2) * gy * gy;
                const float_t mx_hat = mom_x[i] / b1_correction;
                const float_t my_hat = mom_y[i] / b1_correction;
                const float_t vx_hat = second_x[i] / b2_correction;
                const float_t vy_hat = second_y[i] / b2_correction;
                step_x = 0.70f * lr * mx_hat / (std::sqrt(vx_hat) + 1.0e-6f);
                step_y = 0.70f * lr * my_hat / (std::sqrt(vy_hat) + 1.0e-6f);
                const float_t max_step = 4.0f * lr;
                step_x = std::max(-max_step, std::min(max_step, step_x));
                step_y = std::max(-max_step, std::min(max_step, step_y));
            } else {
                if (cfg.use_adaptive_restart &&
                    mom_x[i] * gx + mom_y[i] * gy < 0.0f) {
                    mom_x[i] = 0.0f;
                    mom_y[i] = 0.0f;
                }
                mom_x[i] = momentum * mom_x[i] + (1.0f - momentum) * gx;
                mom_y[i] = momentum * mom_y[i] + (1.0f - momentum) * gy;
                step_x = lr * mom_x[i];
                step_y = lr * mom_y[i];
            }

            cells[i].x -= step_x;
            cells[i].y -= step_y;
            const float_t margin_x = allow_margin ? cells[i].width : 0.0f;
            const float_t margin_y = allow_margin ? cells[i].height : 0.0f;
            cells[i].x = std::max(chip_xl - margin_x,
                                  std::min(chip_xh + margin_x, cells[i].x));
            cells[i].y = std::max(chip_yl - margin_y,
                                  std::min(chip_yh + margin_y, cells[i].y));
        }
    };

    auto compute_density_phase_wirelength = [&]() {
        return compute_hpwl_and_gradient(cells, nets, grad_x, grad_y);
    };

    // ---- CSV convergence log ----
    // Logs: iter, phase, hpwl, objective, avg_overflow_real, dpen, lambda, g_rms, lr
    std::ofstream csv_log("nsp_convergence.csv");
    csv_log << "iter,phase,hpwl,objective,avg_overflow,dpen,lambda,g_rms,lr\n";

    // Helper: compute real avg_overflow from bin density data
    auto compute_real_overflow = [](const BinGrid& grid, const std::vector<float_t>& bin_area_sum,
                                     int nbins) -> float_t {
        float_t bin_area = grid.bin_w * grid.bin_h;
        double excess_area = 0.0;
        for (int b = 0; b < nbins; ++b) {
            float_t rho = bin_area_sum[b] / bin_area;
            float_t over = rho - grid.target_density;
            if (over > 0.0f) excess_area += (double)over * bin_area;
        }
        return grid.available_area > 0.0f
            ? (float_t)(excess_area / grid.available_area) : 0.0f;
    };

    // ================================================================
    // Phase 0: HPWL-only optimization (preserve QP quality)
    // ================================================================
    int n0 = cfg.p0_iters;
    printf("\n--- Phase 0: HPWL-only (%d iters, λ=0) ---\n", n0);
    float_t lambda = 0.0f;
    float_t base_step = (chip_w + chip_h) * 0.005f;
    float_t p0_best_hpwl = std::numeric_limits<float_t>::infinity();
    std::vector<Cell> p0_best_cells;

    for (int iter = 0; iter < n0; ++iter, ++total_iters) {
        hpwl = compute_hpwl_and_gradient(cells, nets, grad_x, grad_y);

        // The QP/ePlace input can already be better than any subgradient
        // iterate. Do not hand a degraded wirelength-only state to spreading.
        if (hpwl < p0_best_hpwl) {
            p0_best_hpwl = hpwl;
            p0_best_cells = cells;
        }

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
    if (!p0_best_cells.empty()) {
        cells = std::move(p0_best_cells);
        hpwl = p0_best_hpwl;
        printf("[Solver] Restored best Phase 0 state: HPWL=%.1f\n", hpwl);
    }

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

    BaselineLambdaConfig baseline_cfg;
    baseline_cfg.interval = cfg.lambda_update_interval;
    baseline_cfg.hpwl_baseline = cfg.hpwl_baseline;
    baseline_cfg.overflow_baseline = cfg.overflow_baseline;
    baseline_cfg.increase_factor = cfg.lambda_increase_factor;
    baseline_cfg.decrease_factor = cfg.lambda_decrease_factor;
    baseline_cfg.hysteresis = cfg.lambda_ratio_hysteresis;
    baseline_cfg.lambda_min = cfg.lambda_control_min;
    baseline_cfg.lambda_max = cfg.lambda_control_max;
    baseline_cfg.guarded = cfg.use_lambda_guard;
    baseline_cfg.recovery_overflow = cfg.lambda_guard_low;
    baseline_cfg.spread_lambda_floor = cfg.lambda_guard_floor;
    baseline_cfg.interval = cfg.lambda_guard_interval;
    BaselineLambdaController guard_controller(baseline_cfg);
    BaselineLambdaConfig early_cfg = baseline_cfg;
    early_cfg.guarded = false;
    early_cfg.interval = cfg.lambda_update_interval;
    BaselineLambdaController early_controller(early_cfg);
    TrajectoryLambdaConfig trajectory_cfg;
    trajectory_cfg.interval = cfg.lambda_guard_interval;
    const int fine_budget = std::max(1, cfg.max_iters - cfg.p0_iters -
        cfg.coarse_iters - cfg.medium_iters);
    trajectory_cfg.horizon_steps = std::min(cfg.lambda_trajectory_horizon,
        std::max(100, fine_budget - 150));
    trajectory_cfg.hpwl_baseline = cfg.hpwl_baseline;
    trajectory_cfg.target_overflow = cfg.lambda_trajectory_target;
    trajectory_cfg.lambda_min = cfg.lambda_control_min;
    trajectory_cfg.lambda_max = cfg.lambda_control_max;
    trajectory_cfg.startup_floor = cfg.lambda_guard_floor;
    TrajectoryLambdaController trajectory_controller(trajectory_cfg);
    std::ofstream lambda_strategy_log("lambda_strategy.csv");
    lambda_strategy_log
        << "global_iteration,phase,density_step,hpwl,overflow,hpwl_ratio,"
           "overflow_ratio,lambda_before,lambda_after,action,controller,"
           "desired_overflow,error,trend_error,integral_error\n";

    auto update_lambda_from_baselines = [&](int phase, float_t current_hpwl,
                                             float_t current_overflow) {
        if (cfg.use_trajectory_lambda && phase == 3) {
            TrajectoryLambdaUpdate update = trajectory_controller.observe(
                lambda, current_hpwl, current_overflow);
            if (!update.evaluated) return;
            lambda = update.lambda_after;
            lambda_strategy_log
                << total_iters + 1 << ',' << phase << ',' << update.density_step << ','
                << current_hpwl << ',' << current_overflow << ','
                << update.hpwl_ratio << ',' << update.overflow_ratio << ','
                << update.lambda_before << ',' << update.lambda_after << ','
                << baseline_lambda_action_name(update.action) << ",trajectory,"
                << update.desired_overflow << ',' << update.error << ','
                << update.trend_error << ',' << update.integral_error << '\n';
            printf("  [LambdaTrajectory step %d] desired=%.4f error=%+.3f "
                   "trend=%+.3f action=%s lambda_control %.6g -> %.6g "
                   "lambda_eff_next=%.4e\n",
                   update.density_step, update.desired_overflow, update.error,
                   update.trend_error, baseline_lambda_action_name(update.action),
                   update.lambda_before, update.lambda_after,
                   update.lambda_after * density_weight_scale);
            return;
        }
        BaselineLambdaController& controller =
            (cfg.use_lambda_guard && phase == 3) ? guard_controller : early_controller;
        BaselineLambdaUpdate update = controller.observe(
            lambda, current_hpwl, current_overflow);
        if (!update.evaluated) return;

        lambda = update.lambda_after;
        lambda_strategy_log
            << total_iters + 1 << ',' << phase << ',' << update.density_step << ','
            << current_hpwl << ',' << current_overflow << ','
            << update.hpwl_ratio << ',' << update.overflow_ratio << ','
            << update.lambda_before << ',' << update.lambda_after << ','
            << baseline_lambda_action_name(update.action)
            << ",baseline,,,,\n";
        printf("  [LambdaBaseline step %d] hpwl_ratio=%.4f overflow_ratio=%.4f "
               "action=%s lambda_control %.6g -> %.6g lambda_eff_next=%.4e\n",
               update.density_step, update.hpwl_ratio, update.overflow_ratio,
               baseline_lambda_action_name(update.action), update.lambda_before,
               update.lambda_after, update.lambda_after * density_weight_scale);
    };

    if (!cfg.use_bayesian_opt) {
        printf("[LambdaBaseline] every %d density steps: HPWL baseline=%.1f, "
               "overflow baseline=%.4f, ratio=baseline/current, hysteresis=%.1f%%\n",
               baseline_cfg.interval, baseline_cfg.hpwl_baseline,
               baseline_cfg.overflow_baseline, baseline_cfg.hysteresis * 100.0f);
        if (cfg.use_trajectory_lambda)
            printf("[LambdaTrajectory] fine target=%.4f horizon=%d steps, "
                   "PI-D trajectory feedback enabled\n",
                   trajectory_cfg.target_overflow, trajectory_cfg.horizon_steps);
    }

    printf("\n--- Phase 1: Coarse Spread (%d iters) ---\n", n1);
    reset_density_optimizer();
    for (int iter = 0; iter < n1; ++iter, ++total_iters) {
        const BinGrid& grid = grid_coarse;
        bin_buf.resize(grid.nx * grid.ny);

        hpwl = compute_density_phase_wirelength();
        std::fill(density_grad_x.begin(), density_grad_x.end(), 0.0f);
        std::fill(density_grad_y.begin(), density_grad_y.end(), 0.0f);
        dpen = density_compute_and_gradient(
            grid, cells, density_grad_x, density_grad_y, bin_buf);
        const float_t effective_lambda = combine_density_gradient(lambda);

        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            double gn = (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
            gsum2 += gn;
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        float_t progress = (float_t)iter / n1;
        float_t lr = base_step / std::sqrt(progress * 4.0f + 1.0f);

        apply_density_step(g_rms, lr, cfg.momentum, true);

        if (hpwl < hpwl_best) hpwl_best = hpwl;

        float_t real_overflow = compute_real_overflow(grid, bin_buf, grid.nx * grid.ny);
        float_t avg_over = real_overflow;

        if (cfg.use_bayesian_opt) {
            lambda = (float_t)lambda_bo.step((double)lambda, (double)hpwl,
                                              (double)real_overflow, (double)target_overflow);
        } else {
            update_lambda_from_baselines(1, hpwl, real_overflow);
        }

        float_t objective = hpwl + effective_lambda * dpen;

        csv_log << iter << ",1," << hpwl << "," << objective << ","
                << real_overflow << "," << dpen << "," << effective_lambda << ","
                << g_rms << "," << lr << "\n";
        visualization_snapshot(cfg.visualization, total_iters + 1, "coarse", hpwl,
                               real_overflow, dpen, effective_lambda, cells);

        if (iter % 30 == 0 || iter == n1 - 1)
            printf("  [c%3d] HPWL=%.1f Dpen=%.1f lambda_eff=%.4e lr=%.2f g_rms=%.3f ov=%.4f\n",
                    iter, hpwl, dpen, effective_lambda, lr, g_rms, avg_over);
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
    density_weight_scale = -1.0f;

    printf("\n--- Phase 2: Medium Balance (%d iters) ---\n", n2);
    reset_density_optimizer();
    for (int iter = 0; iter < n2; ++iter, ++total_iters) {
        const BinGrid& grid = grid_medium;
        bin_buf.resize(grid.nx * grid.ny);

        hpwl = compute_density_phase_wirelength();
        std::fill(density_grad_x.begin(), density_grad_x.end(), 0.0f);
        std::fill(density_grad_y.begin(), density_grad_y.end(), 0.0f);
        dpen = density_compute_and_gradient(
            grid, cells, density_grad_x, density_grad_y, bin_buf);
        const float_t effective_lambda = combine_density_gradient(lambda);

        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            double gn = (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
            gsum2 += gn;
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        float_t progress = (float_t)iter / n2;
        float_t lr = base_step / std::sqrt(progress * 4.0f + 1.0f);

        apply_density_step(g_rms, lr, 0.80f, true);

        if (hpwl < hpwl_best) hpwl_best = hpwl;

        float_t real_overflow = compute_real_overflow(grid, bin_buf, grid.nx * grid.ny);
        float_t avg_over = real_overflow;

        if (cfg.use_bayesian_opt) {
            lambda = (float_t)lambda_bo.step((double)lambda, (double)hpwl,
                                              (double)real_overflow, (double)target_overflow);
        } else {
            update_lambda_from_baselines(2, hpwl, real_overflow);
        }

        float_t objective = hpwl + effective_lambda * dpen;

        csv_log << iter << ",2," << hpwl << "," << objective << ","
                << real_overflow << "," << dpen << "," << effective_lambda << ","
                << g_rms << "," << lr << "\n";
        visualization_snapshot(cfg.visualization, total_iters + 1, "medium", hpwl,
                               real_overflow, dpen, effective_lambda, cells);

        if (iter % 50 == 0 || iter == n2 - 1)
            printf("  [m%3d] HPWL=%.1f Dpen=%.1f lambda_eff=%.4e lr=%.2f g_rms=%.3f ov=%.4f\n",
                    iter, hpwl, dpen, effective_lambda, lr, g_rms, avg_over);
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
    density_weight_scale = -1.0f;
    float_t best_feasible_hpwl = std::numeric_limits<float_t>::infinity();
    float_t best_feasible_overflow = 0.0f;
    int best_feasible_iter = -1;
    std::vector<Cell> best_feasible_cells;

    printf("\n--- Phase 3: Fine Refinement (%d iters) ---\n", n3);
    printf("[Feasibility] best-state overflow limit=%.6f\n", cfg.overflow_limit);
    if (cfg.use_projected_hpwl)
        printf("[ProjectedHPWL] active in overflow band [%.6f, %.6f]\n",
               std::max(0.0f, cfg.overflow_baseline - cfg.projection_overflow_margin),
               cfg.overflow_baseline + cfg.projection_overflow_margin);
    if (cfg.use_adaptive_restart)
        printf("[AdaptiveRestart] local heavy-ball conflict restart enabled\n");
    if (cfg.use_hpwl_bundle)
        printf("[HPWLBundle] cuts=%d interval=%d prox=%.3f current_mix=%.3f%s%s\n",
               cfg.hpwl_bundle_size, cfg.hpwl_bundle_interval,
               cfg.hpwl_bundle_prox_scale, cfg.hpwl_bundle_current_mix,
               cfg.hpwl_bundle_adaptive_cuts ? " adaptive-cuts" : "",
               cfg.hpwl_bundle_extrema_mix ? " extrema-mix" : "");
    reset_density_optimizer();
    for (int iter = 0; iter < n3; ++iter, ++total_iters) {
        const BinGrid& grid = grid_fine;
        bin_buf.resize(grid.nx * grid.ny);

        hpwl = compute_density_phase_wirelength();
        std::fill(density_grad_x.begin(), density_grad_x.end(), 0.0f);
        std::fill(density_grad_y.begin(), density_grad_y.end(), 0.0f);
        dpen = density_compute_and_gradient(
            grid, cells, density_grad_x, density_grad_y, bin_buf);
        float_t real_overflow = compute_real_overflow(grid, bin_buf, grid.nx * grid.ny);
        if (cfg.use_projected_hpwl &&
            real_overflow >= cfg.overflow_baseline - cfg.projection_overflow_margin &&
            real_overflow <= cfg.overflow_baseline + cfg.projection_overflow_margin)
            project_hpwl_from_density_ascent();

        float_t progress = (float_t)iter / n3;
        float_t lr = base_step / std::sqrt(progress * 8.0f + 1.0f);
        const int cooldown_start = n3 > 800
            ? std::min(1800, n3 - 800) : n3;
        if (iter > cooldown_start) {
            const float_t cooldown_progress = (float_t)(iter - cooldown_start) /
                std::max(1, n3 - 1 - cooldown_start);
            const float_t cooldown_floor = cfg.fine_cooldown_floor;
            lr *= std::max(cooldown_floor,
                           1.0f - (1.0f - cooldown_floor) * cooldown_progress);
        }

        double hpwl_gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            const int i = movable_ids[idx];
            hpwl_gsum2 += (double)grad_x[i] * grad_x[i] +
                          (double)grad_y[i] * grad_y[i];
        }
        const float_t hpwl_g_rms =
            (float_t)std::sqrt(hpwl_gsum2 / Nm + RMS_EPS);
        HpwlBundleStats bundle_stats =
            apply_hpwl_bundle(iter, hpwl, lr, hpwl_g_rms);
        const float_t effective_lambda = combine_density_gradient(lambda);

        double gsum2 = 0.0;
        for (int idx = 0; idx < Nm; ++idx) {
            int i = movable_ids[idx];
            double gn = (double)grad_x[i]*grad_x[i] + (double)grad_y[i]*grad_y[i];
            gsum2 += gn;
        }
        float_t g_rms = (float_t)std::sqrt(gsum2 / Nm + RMS_EPS);

        if (hpwl < hpwl_best) hpwl_best = hpwl;

        float_t avg_over = real_overflow;

        if (real_overflow <= cfg.overflow_limit && hpwl < best_feasible_hpwl) {
            best_feasible_hpwl = hpwl;
            best_feasible_overflow = real_overflow;
            best_feasible_iter = iter;
            best_feasible_cells = cells;
        }

        if (cfg.use_bayesian_opt) {
            lambda = (float_t)lambda_bo.step((double)lambda, (double)hpwl,
                                              (double)real_overflow, (double)target_overflow);
        } else {
            update_lambda_from_baselines(3, hpwl, real_overflow);
        }

        float_t objective = hpwl + effective_lambda * dpen;

        csv_log << iter << ",3," << hpwl << "," << objective << ","
                << real_overflow << "," << dpen << "," << effective_lambda << ","
                << g_rms << "," << lr << "\n";
        visualization_snapshot(cfg.visualization, total_iters + 1, "fine", hpwl,
                               real_overflow, dpen, effective_lambda, cells);
        if (bundle_stats.updated) {
            bundle_log << total_iters + 1 << ',' << iter << ',' << bundle_stats.cuts << ','
                       << bundle_stats.tau << ',' << bundle_stats.newest_weight << ','
                       << bundle_stats.model_error << ',' << bundle_stats.norm_ratio << ','
                       << bundle_stats.entropy << ',' << bundle_stats.effective_current_mix << ','
                       << bundle_stats.trigger_cosine << ',' << bundle_stats.cut_age << '\n';
        }

        if (iter % 100 == 0 || iter == n3 - 1)
            printf("  [f%3d] HPWL=%.1f Dpen=%.1f lambda_eff=%.4e lr=%.2f g_rms=%.3f ov=%.4f\n",
                    iter, hpwl, dpen, effective_lambda, lr, g_rms, avg_over);

        apply_density_step(g_rms, lr, cfg.fine_momentum, false);
    }

    csv_log.close();
    lambda_strategy_log.close();
    if (bundle_log.is_open()) bundle_log.close();
    printf("[Solver] Convergence log saved to nsp_convergence.csv\n");

    if (!best_feasible_cells.empty()) {
        cells = std::move(best_feasible_cells);
        printf("[Solver] Restored best feasible fine-grid state: iter=%d "
               "HPWL=%.1f overflow=%.6f\n",
               best_feasible_iter, best_feasible_hpwl, best_feasible_overflow);
    }

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
