#include "nonsmooth.h"

#include "bookshelf.h"
#include "legalizer.h"
#include "obstacle.h"
#include "wirelength.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace dpcpp {
namespace {

struct BundleCut {
    Real hpwl = 0.0;
    std::vector<Real> x, y, gx, gy;
};

struct BundleStats {
    int cuts = 0;
    Real newest_weight = 1.0;
    Real model_error = 0.0;
    Real cosine = 1.0;
};

std::vector<Real> capture(const Database& db, const std::vector<Filler>& fillers) {
    const std::size_t n = db.movable_ids.size() + fillers.size();
    const std::size_t movable_count = db.movable_ids.size();
    std::vector<Real> values(2 * n);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        if (static_cast<std::size_t>(i) < movable_count) {
            const Node& node = db.nodes[db.movable_ids[i]];
            values[i] = node.x;
            values[n + i] = node.y;
        } else {
            const std::size_t filler_index = i - movable_count;
            values[i] = fillers[filler_index].x;
            values[n + i] = fillers[filler_index].y;
        }
    }
    return values;
}

void capture_into(const Database& db, const std::vector<Filler>& fillers,
                  std::vector<Real>& values) {
    const std::size_t n = db.movable_ids.size() + fillers.size();
    const std::size_t movable_count = db.movable_ids.size();
    values.resize(2 * n);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        if (static_cast<std::size_t>(i) < movable_count) {
            const Node& node = db.nodes[db.movable_ids[i]];
            values[i] = node.x;
            values[n + i] = node.y;
        } else {
            const std::size_t filler_index = i - movable_count;
            values[i] = fillers[filler_index].x;
            values[n + i] = fillers[filler_index].y;
        }
    }
}

void apply(Database& db, std::vector<Filler>& fillers,
           const std::vector<Real>& values) {
    const std::size_t n = db.movable_ids.size() + fillers.size();
    const std::size_t movable_count = db.movable_ids.size();
    if (values.size() != 2 * n) throw std::runtime_error("position vector size mismatch");
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        if (static_cast<std::size_t>(i) < movable_count) {
            Node& node = db.nodes[db.movable_ids[i]];
            node.x = std::clamp(values[i], db.xl + 0.5 * node.width,
                                db.xh - 0.5 * node.width);
            node.y = std::clamp(values[n + i], db.yl + 0.5 * node.height,
                                db.yh - 0.5 * node.height);
        } else {
            const std::size_t filler_index = i - movable_count;
            Filler& filler = fillers[filler_index];
            filler.x = std::clamp(values[i], db.xl + 0.5 * filler.width,
                                  db.xh - 0.5 * filler.width);
            filler.y = std::clamp(values[n + i], db.yl + 0.5 * filler.height,
                                  db.yh - 0.5 * filler.height);
        }
    }
}

Real l1_norm(const std::vector<Real>& values) {
    Real sum = 0.0;
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (int i = 0; i < static_cast<int>(values.size()); ++i)
        sum += std::abs(values[i]);
    return sum;
}

Real rms_movable(const Database& db, const std::vector<Real>& gx,
                 const std::vector<Real>& gy) {
    Real sum = 0.0;
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (int i = 0; i < static_cast<int>(db.movable_ids.size()); ++i) {
        const int id = db.movable_ids[i];
        sum += gx[id] * gx[id] + gy[id] * gy[id];
    }
    return std::sqrt(sum / std::max<std::size_t>(1, 2 * db.movable_ids.size()));
}

void add_gp_noise(Database& db, std::vector<Filler>& fillers,
                  std::uint64_t seed, Real ratio) {
    if (ratio <= 0.0) return;
    std::mt19937_64 generator(seed ^ 0xd1b54a32d192ed03ULL);
    std::uniform_real_distribution<Real> unit(-0.5, 0.5);
    for (int id : db.movable_ids) {
        Node& node = db.nodes[id];
        node.x += unit(generator) * node.width * ratio;
        node.y += unit(generator) * node.height * ratio;
    }
    for (Filler& filler : fillers) {
        filler.x += unit(generator) * filler.width * ratio;
        filler.y += unit(generator) * filler.height * ratio;
    }
    clamp_to_region(db, fillers);
}

int nearest_compatible_row(const Database& db, const Node& node) {
    int best = 0;
    Real distance = std::numeric_limits<Real>::infinity();
    const Real bottom = node.y - 0.5 * node.height;
    for (int r = 0; r < static_cast<int>(db.rows.size()); ++r) {
        if (node.height > db.rows[r].height + 1.0e-6) continue;
        const Real current = std::abs(bottom - db.rows[r].y);
        if (current < distance) {
            distance = current;
            best = r;
        }
    }
    return best;
}

void add_row_proximity_force(const Database& db, Real strength,
                             std::vector<Real>& gy) {
    if (strength <= 0.0 || db.rows.empty()) return;
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Row& row = db.rows[nearest_compatible_row(db, node)];
        const Real target = row.y + 0.5 * node.height;
        const Real normalized = (node.y - target) /
            std::max<Real>(row.height, 1.0e-9);
        gy[id] += strength * std::clamp(normalized, -1.0, 1.0) *
                  std::max(1, db.node_pin_weight[id]);
    }
}

void sample_exact_subgradient(Database& db, int samples, Real radius,
                              std::uint64_t seed, int degree_limit,
                              std::vector<Real>& gx, std::vector<Real>& gy,
                              int& evaluations) {
    if (samples <= 0 || radius <= 0.0) return;
    std::vector<std::pair<Real, Real>> original(db.nodes.size());
    for (int id : db.movable_ids)
        original[id] = {db.nodes[id].x, db.nodes[id].y};
    std::vector<Real> aggregate_x = gx;
    std::vector<Real> aggregate_y = gy;
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<Real> perturb(-radius, radius);
    for (int sample = 0; sample < samples; ++sample) {
        for (int id : db.movable_ids) {
            Node& node = db.nodes[id];
            node.x = std::clamp(original[id].first + perturb(generator),
                                db.xl + 0.5 * node.width,
                                db.xh - 0.5 * node.width);
            node.y = std::clamp(original[id].second + perturb(generator),
                                db.yl + 0.5 * node.height,
                                db.yh - 0.5 * node.height);
        }
        std::vector<Real> sample_x, sample_y;
        exact_hpwl_subgradient(db, degree_limit, &sample_x, &sample_y);
        ++evaluations;
        for (int id : db.movable_ids) {
            aggregate_x[id] += sample_x[id];
            aggregate_y[id] += sample_y[id];
        }
    }
    for (int id : db.movable_ids) {
        db.nodes[id].x = original[id].first;
        db.nodes[id].y = original[id].second;
        gx[id] = aggregate_x[id] / static_cast<Real>(samples + 1);
        gy[id] = aggregate_y[id] / static_cast<Real>(samples + 1);
    }
}

bool quick_legalize(const Database& db, int detailed_passes,
                    Database& legalized, Real& legal_hpwl) {
    legalized = db;
    LegalizeConfig config;
    config.detailed_passes = std::max(0, detailed_passes);
    config.run_dreamplace_detailed = false;
    config.detailed_outer_rounds = 1;
    try {
        const LegalizeResult result = legalize_and_refine(legalized, config);
        legal_hpwl = result.hpwl_after_detailed;
        return result.legality.legal;
    } catch (const std::exception& error) {
        std::cout << "[QuickLegal] failed=" << error.what() << '\n';
        return false;
    }
}

int initial_grid_size(int target, int requested_minimum) {
    if (target <= requested_minimum) return target;
    int value = 1;
    while (value * 2 <= target && value * 2 <= requested_minimum) value *= 2;
    return value;
}

Real vector_cosine(const Database& db, const std::vector<Real>& ax,
                   const std::vector<Real>& ay, const BundleCut& b) {
    Real dot = 0.0, an = 0.0, bn = 0.0;
    for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
        const int id = db.movable_ids[i];
        dot += ax[id] * b.gx[i] + ay[id] * b.gy[i];
        an += ax[id] * ax[id] + ay[id] * ay[id];
        bn += b.gx[i] * b.gx[i] + b.gy[i] * b.gy[i];
    }
    return dot / std::max<Real>(std::sqrt(an * bn), 1.0e-30);
}

void append_bundle_cut(const Database& db, Real hpwl,
                       const std::vector<Real>& gx,
                       const std::vector<Real>& gy,
                       const GlobalPlaceConfig& config,
                       std::vector<BundleCut>& cuts) {
    BundleCut cut;
    cut.hpwl = hpwl;
    const std::size_t n = db.movable_ids.size();
    cut.x.resize(n);
    cut.y.resize(n);
    cut.gx.resize(n);
    cut.gy.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const int id = db.movable_ids[i];
        cut.x[i] = db.nodes[id].x;
        cut.y[i] = db.nodes[id].y;
        cut.gx[i] = gx[id];
        cut.gy[i] = gy[id];
    }
    if (static_cast<int>(cuts.size()) >= config.bundle_size) cuts.erase(cuts.begin());
    cuts.push_back(std::move(cut));
}

BundleStats apply_bundle(const Database& db, Real hpwl, Real learning_rate,
                         int iteration, const GlobalPlaceConfig& config,
                         std::vector<BundleCut>& cuts, int& last_cut,
                         Real prox_multiplier,
                         std::vector<Real>& gx, std::vector<Real>& gy) {
    BundleStats stats;
    Real cosine = cuts.empty() ? 1.0 : vector_cosine(db, gx, gy, cuts.back());
    const int age = iteration - last_cut;
    bool insert = cuts.empty() || age >= config.bundle_interval;
    if (!cuts.empty() && age >= config.bundle_min_interval &&
        cosine <= config.bundle_cosine_trigger) insert = true;
    if (insert) {
        append_bundle_cut(db, hpwl, gx, gy, config, cuts);
        last_cut = iteration;
    }

    const int kcount = static_cast<int>(cuts.size());
    stats.cuts = kcount;
    stats.cosine = cosine;
    if (kcount == 0) return stats;
    std::vector<Real> beta(kcount, 0.0), gram(kcount * kcount, 0.0);
    for (int j = 0; j < kcount; ++j) {
        Real shift = 0.0;
        for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
            const int id = db.movable_ids[i];
            shift += cuts[j].gx[i] * (db.nodes[id].x - cuts[j].x[i]) +
                     cuts[j].gy[i] * (db.nodes[id].y - cuts[j].y[i]);
        }
        beta[j] = cuts[j].hpwl + shift;
        for (int k = 0; k <= j; ++k) {
            Real value = 0.0;
            for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
                value += cuts[j].gx[i] * cuts[k].gx[i] +
                         cuts[j].gy[i] * cuts[k].gy[i];
            }
            gram[j * kcount + k] = gram[k * kcount + j] = value;
        }
    }
    const Real tau = config.bundle_prox_scale * prox_multiplier * learning_rate /
                     std::max<Real>(rms_movable(db, gx, gy), 1.0e-9);
    std::vector<Real> alpha(kcount, 0.0);
    alpha.back() = 1.0;
    for (int fw = 0; fw < 20; ++fw) {
        std::vector<Real> dual_gradient(kcount, 0.0);
        int best = 0;
        for (int j = 0; j < kcount; ++j) {
            Real q = 0.0;
            for (int k = 0; k < kcount; ++k) q += gram[j * kcount + k] * alpha[k];
            dual_gradient[j] = beta[j] - tau * q;
            if (dual_gradient[j] > dual_gradient[best]) best = j;
        }
        std::vector<Real> direction(kcount);
        for (int j = 0; j < kcount; ++j) direction[j] = -alpha[j];
        direction[best] += 1.0;
        Real numerator = 0.0, curvature = 0.0;
        for (int j = 0; j < kcount; ++j) {
            numerator += dual_gradient[j] * direction[j];
            for (int k = 0; k < kcount; ++k)
                curvature += direction[j] * gram[j * kcount + k] * direction[k];
        }
        if (numerator <= 1.0e-7) break;
        const Real fraction = curvature > 1.0e-20
            ? std::min<Real>(1.0, numerator / (tau * curvature)) : 1.0;
        for (int j = 0; j < kcount; ++j) alpha[j] += fraction * direction[j];
    }

    std::vector<Real> aggregate_x(db.nodes.size(), 0.0);
    std::vector<Real> aggregate_y(db.nodes.size(), 0.0);
    Real aggregate_beta = 0.0;
    for (int j = 0; j < kcount; ++j) {
        aggregate_beta += alpha[j] * beta[j];
        for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
            const int id = db.movable_ids[i];
            aggregate_x[id] += alpha[j] * cuts[j].gx[i];
            aggregate_y[id] += alpha[j] * cuts[j].gy[i];
        }
    }
    const Real mix = std::clamp(config.bundle_current_mix, 0.0, 1.0);
    for (int id : db.movable_ids) {
        gx[id] = mix * gx[id] + (1.0 - mix) * aggregate_x[id];
        gy[id] = mix * gy[id] + (1.0 - mix) * aggregate_y[id];
    }
    stats.newest_weight = alpha.back();
    stats.model_error = std::max<Real>(0.0, hpwl - aggregate_beta);
    return stats;
}

BundleStats apply_grouped_bundle(
    const Database& db, const std::vector<Real>& group_hpwl, Real learning_rate,
    int iteration, const GlobalPlaceConfig& config,
    std::vector<std::vector<BundleCut>>& grouped_cuts,
    std::vector<int>& grouped_last_cut, Real prox_multiplier,
    std::vector<std::vector<Real>>& group_gx,
    std::vector<std::vector<Real>>& group_gy,
    std::vector<Real>& aggregate_gx, std::vector<Real>& aggregate_gy) {
    const int groups = static_cast<int>(group_hpwl.size());
    BundleStats total;
    total.newest_weight = 0.0;
    total.cosine = 0.0;
    GlobalPlaceConfig group_config = config;
    group_config.bundle_size = config.bundle_group_size;
    for (int group = 0; group < groups; ++group) {
        const BundleStats current = apply_bundle(
            db, group_hpwl[group], learning_rate, iteration, group_config,
            grouped_cuts[group], grouped_last_cut[group], prox_multiplier,
            group_gx[group], group_gy[group]);
        total.cuts += current.cuts;
        total.newest_weight += current.newest_weight;
        total.model_error += current.model_error;
        total.cosine += current.cosine;
    }
    total.newest_weight /= std::max(1, groups);
    total.cosine /= std::max(1, groups);
    aggregate_gx.assign(db.nodes.size(), 0.0);
    aggregate_gy.assign(db.nodes.size(), 0.0);
    for (int group = 0; group < groups; ++group) {
        for (int id : db.movable_ids) {
            aggregate_gx[id] += group_gx[group][id];
            aggregate_gy[id] += group_gy[group][id];
        }
    }
    return total;
}

Real update_lambda_control(Real eta, LambdaPolicy policy, int density_step,
                           Real hpwl, Real overflow, Real previous_hpwl,
                           Real initial_overflow, Real& previous_overflow,
                           Real& integral_error,
                           const GlobalPlaceConfig& config) {
    if (policy == LambdaPolicy::Dreamplace) {
        const Real delta = hpwl - previous_hpwl;
        Real multiplier = 1.0;
        if (delta < 0.0) {
            multiplier = 1.05 * std::max(std::pow(0.9999, density_step), 0.98);
        } else {
            multiplier = std::clamp(
                std::pow(1.05, 1.0 - delta / config.reference_hpwl_delta),
                0.95, 1.05);
        }
        return std::clamp(eta * multiplier,
                          config.lambda_control_min, config.lambda_control_max);
    }
    if (config.lambda_update_interval <= 0 ||
        density_step % config.lambda_update_interval != 0) return eta;
    if (policy == LambdaPolicy::Ratio) {
        const Real hpwl_target = config.hpwl_baseline > 0.0
            ? config.hpwl_baseline : hpwl;
        const Real hpwl_ratio = hpwl / std::max<Real>(hpwl_target, 1.0);
        const Real overflow_ratio = overflow /
                                    std::max<Real>(config.overflow_baseline, 1.0e-9);
        Real factor = 1.0;
        // Lambda multiplies density: excess overflow must increase it, while
        // an HPWL-dominated feasible state can release density pressure.
        if (overflow_ratio > 1.05 * hpwl_ratio) factor = 1.5;
        else if (hpwl_ratio > 1.05 * overflow_ratio) factor = 2.0 / 3.0;
        if (overflow > 1.05 * config.stop_overflow) factor = std::max(factor, 1.25);
        return std::clamp(eta * factor,
                          config.lambda_control_min, config.lambda_control_max);
    }

    const Real progress = std::min<Real>(1.0, density_step /
        static_cast<Real>(std::max(1, config.lambda_trajectory_horizon)));
    const Real smoothstep = progress * progress * (3.0 - 2.0 * progress);
    const Real desired = config.stop_overflow +
        (initial_overflow - config.stop_overflow) * (1.0 - smoothstep);
    const Real scale = std::max<Real>(config.stop_overflow, 1.0e-4);
    const Real error = std::clamp((overflow - desired) / scale, -2.0, 4.0);
    const Real desired_drop = (previous_overflow - desired) / scale;
    const Real actual_drop = (previous_overflow - overflow) / scale;
    const Real trend_error = std::clamp(desired_drop - actual_drop, -2.0, 2.0);
    if (progress < 1.0) integral_error = std::clamp(integral_error + error, -3.0, 3.0);
    else integral_error = 0.0;
    // The feedback is evaluated only once per interval.  Express the gains
    // as an equivalent per-iteration log change and accumulate them across
    // the interval; otherwise a 25-step controller ramps about 25x too
    // slowly compared with DREAMPlace's per-iteration lambda update.
    Real log_delta_per_iteration = 0.045 + 0.08 * error +
        0.003 * integral_error + 0.04 * trend_error;
    if (overflow <= config.stop_overflow && config.hpwl_baseline > 0.0) {
        const Real pressure = std::max<Real>(0.0, hpwl / config.hpwl_baseline - 1.0);
        log_delta_per_iteration -= 0.01 * pressure;
    }
    log_delta_per_iteration = std::clamp(
        log_delta_per_iteration, -0.05, 0.06);
    const Real log_delta = std::clamp(
        log_delta_per_iteration * config.lambda_update_interval, -1.0, 1.25);
    previous_overflow = overflow;
    return std::clamp(eta * std::exp(log_delta),
                      config.lambda_control_min, config.lambda_control_max);
}

Real coordinate_step(GlobalOptimizer optimizer, Real gradient, Real learning_rate,
                     const GlobalPlaceConfig& config, Real& first,
                     Real& second, Real& maximum_second,
                     Real first_bias_denominator,
                     Real second_bias_denominator) {
    if (optimizer == GlobalOptimizer::HeavyBall) {
        first = config.momentum * first + (1.0 - config.momentum) * gradient;
        return learning_rate * first;
    }
    if (optimizer == GlobalOptimizer::AdaGrad) {
        second += gradient * gradient;
        return std::clamp(learning_rate * gradient /
            (std::sqrt(second) + config.epsilon),
            -config.max_step_multiplier * learning_rate,
             config.max_step_multiplier * learning_rate);
    }
    first = config.beta1 * first + (1.0 - config.beta1) * gradient;
    second = config.beta2 * second + (1.0 - config.beta2) * gradient * gradient;
    const Real first_hat = first / first_bias_denominator;
    const Real second_hat = second / second_bias_denominator;
    Real denominator = second_hat;
    if (optimizer == GlobalOptimizer::AMSGrad) {
        maximum_second = std::max(maximum_second, second_hat);
        denominator = maximum_second;
    }
    return std::clamp(learning_rate * first_hat /
        (std::sqrt(denominator) + config.epsilon),
        -config.max_step_multiplier * learning_rate,
         config.max_step_multiplier * learning_rate);
}

Real feasible_band_lambda(Real eta, Real overflow,
                          const GlobalPlaceConfig& config) {
    const Real upper = config.stop_overflow;
    const Real lower = config.refinement_lower_overflow;
    const Real band = std::max<Real>(upper - lower, 1.0e-4);
    Real log_delta = 0.0;
    if (overflow > upper) {
        log_delta = config.refinement_lambda_gain * (overflow - upper) / band;
    } else if (overflow < lower) {
        log_delta = -config.refinement_lambda_gain * (lower - overflow) / band;
    } else {
        // Release density pressure gradually inside the feasible band.  The
        // release vanishes at its upper edge, so an infeasible step reverses
        // the controller immediately instead of allowing lambda to run away.
        log_delta = -0.05 * config.refinement_lambda_gain *
                    (upper - overflow) / band;
    }
    log_delta = std::clamp(log_delta, -0.05, 0.05);
    return std::clamp(eta * std::exp(log_delta),
                      config.lambda_control_min, config.lambda_control_max);
}

Real progressive_band_dual_lambda(Real eta, int density_step, Real hpwl,
                                  Real overflow, Real initial_overflow,
                                  Real& integral_error,
                                  const GlobalPlaceConfig& config) {
    const Real lower = config.progressive_overflow_lower;
    const Real upper = config.progressive_overflow_upper;
    const Real midpoint = 0.5 * (lower + upper);
    const Real progress = std::min<Real>(
        1.0, density_step /
        static_cast<Real>(std::max(1, config.progressive_density_iterations)));
    const Real smoothstep = progress * progress * (3.0 - 2.0 * progress);
    const Real desired = midpoint + (initial_overflow - midpoint) *
        (1.0 - smoothstep);
    const Real scale = std::max<Real>(upper, 1.0e-4);
    Real error = std::clamp((overflow - desired) / scale, -3.0, 5.0);
    if (progress >= 1.0) {
        if (overflow > upper) error = (overflow - upper) / scale;
        else if (overflow < lower) error = (overflow - lower) / scale;
        else error = 0.0;
    }
    integral_error = std::clamp(integral_error + error, -8.0, 8.0);
    const Real hpwl_excess = std::max<Real>(
        0.0, hpwl / std::max<Real>(config.progressive_hpwl_target, 1.0) - 1.0);
    Real log_delta = config.progressive_dual_kp * error +
        config.progressive_dual_ki * integral_error -
        config.progressive_hpwl_gain * hpwl_excess;
    // Once the target band is reached, lambda is a dual variable: it stays
    // put unless a constraint error or an HPWL budget violation moves it.
    log_delta = std::clamp(log_delta, -0.08, 0.12);
    return std::clamp(eta * std::exp(log_delta),
                      config.lambda_control_min, config.lambda_control_max);
}

void save_global_snapshot(const Database& db, const GlobalPlaceConfig& config,
                          int iteration) {
    if (config.snapshot_every <= 0 || config.snapshot_dir.empty() ||
        iteration % config.snapshot_every != 0) return;
    std::ostringstream name;
    name << "global_" << std::setw(4) << std::setfill('0') << iteration << ".pl";
    write_bookshelf_pl(db, (std::filesystem::path(config.snapshot_dir) /
                            name.str()).string());
}

}  // namespace

GlobalPlaceResult nonsmooth_global_place(Database& db,
                                         std::vector<Filler>& fillers,
                                         const GlobalPlaceConfig& config) {
    const auto start = std::chrono::steady_clock::now();
    int grid_x = config.multilevel_density
        ? initial_grid_size(config.bins_x, config.multilevel_min_bins)
        : config.bins_x;
    int grid_y = config.multilevel_density
        ? initial_grid_size(config.bins_y, config.multilevel_min_bins)
        : config.bins_y;
    int total_grid_transitions = 0;
    for (int x = grid_x, y = grid_y;
         x < config.bins_x || y < config.bins_y;
         x = std::min(config.bins_x, x * 2),
         y = std::min(config.bins_y, y * 2)) {
        ++total_grid_transitions;
    }
    int completed_grid_transitions = 0;
    int last_grid_transition = -1000000;
    const ElectricFieldModel field_model = config.mixed_spectral_field
        ? ElectricFieldModel::MixedSpectral : ElectricFieldModel::FiniteDifference;
    std::unique_ptr<ElectricDensity> density = std::make_unique<ElectricDensity>(
        db, grid_x, grid_y, config.target_density, field_model);
    std::unique_ptr<FixedObstacleField> obstacle_field;
    if (config.progressive_legalization) {
        obstacle_field = std::make_unique<FixedObstacleField>(db);
    }
    add_gp_noise(db, fillers, config.seed, config.gp_noise_ratio);
    std::ofstream metrics;
    if (!config.metrics_path.empty()) {
        std::filesystem::create_directories(
            std::filesystem::path(config.metrics_path).parent_path());
        metrics.open(config.metrics_path);
        metrics << "iteration,exact_hpwl,overflow,max_density,density_energy,"
                   "lambda_base,lambda_control,lambda_effective,learning_rate,"
                   "bundle_cuts,bundle_newest_weight,bundle_model_error,bundle_cosine,"
                   "refinement,optimizer,grid_x,grid_y,row_distance,segment_overflow,"
                   "quick_legal_hpwl,serious_step,active_radius_scale,active_radius,"
                   "adaptive_signal,adaptive_phase";
        if (config.progressive_legalization) {
            metrics << ",stage,obstacle_energy,obstacle_weight,macro_overlap_area,"
                       "macro_overlap_ratio,macro_overlap_cells,"
                       "previous_filter_fraction";
        }
        metrics << '\n';
        metrics << std::setprecision(12);
    }

    const std::size_t n = db.movable_ids.size() + fillers.size();
    std::vector<Real> first(2 * n, 0.0), second(2 * n, 0.0), maximum_second(2 * n, 0.0);
    std::vector<Real> dgx, dgy, fgx, fgy;
    std::vector<Real> wgx, wgy;
    std::vector<Real> gradient(2 * n, 0.0);
    std::vector<Real> previous_positions(2 * n), positions(2 * n);
    std::vector<Real> wire_gradient_refinement(2 * n, 0.0);
    std::vector<Real> density_normal_refinement(2 * n, 0.0);
    std::vector<BundleCut> cuts;
    std::vector<std::vector<BundleCut>> grouped_cuts(
        config.bundle_groups > 1 ? config.bundle_groups : 0);
    std::vector<int> grouped_last_cut(
        config.bundle_groups > 1 ? config.bundle_groups : 0, -1000000);
    int last_cut = -1000000;
    auto clear_bundle_models = [&]() {
        cuts.clear();
        last_cut = -1000000;
        for (auto& group : grouped_cuts) group.clear();
        std::fill(grouped_last_cut.begin(), grouped_last_cut.end(), -1000000);
    };
    Real lambda_base = 0.0;
    Real lambda_control = 1.0;
    Real pending_effective_lambda = 0.0;
    Real previous_hpwl = exact_hpwl(db);
    Real initial_density_overflow = -1.0;
    Real previous_overflow = 0.0;
    Real integral_error = 0.0;
    int density_step = 0;
    bool refinement_active = false;
    int refinement_start = -1;
    Real bundle_prox_multiplier = 1.0;
    Real serious_trial_scale = 1.0;
    Real tangent_trial_scale = 1.0;
    Real active_radius_scale = 1.0;
    Real adaptive_signal_ema = 0.0;
    int adaptive_phase = 0;
    std::deque<Real> adaptive_hpwl_history;
    std::deque<Real> adaptive_overflow_history;
    int last_serious_step = 1;
    Real obstacle_weight_base = 0.0;
    int previous_stage = 0;
    Real last_filter_fraction = 1.0;
    GlobalPlaceResult result;
    double profile_density = 0.0, profile_wire = 0.0;
    double profile_update = 0.0;
    std::vector<Real> best_feasible, best_overflow = capture(db, fillers);
    std::vector<Real> best_legal_raw;
    Metrics best_legal_metrics;
    Real best_legal_hpwl = std::numeric_limits<Real>::infinity();
    Metrics best_overflow_metrics;
    best_overflow_metrics.overflow = std::numeric_limits<Real>::infinity();
    std::vector<Real> best_progressive;
    Metrics best_progressive_metrics;
    Real best_progressive_score = std::numeric_limits<Real>::infinity();

    for (int iteration = 0; iteration < config.iterations; ++iteration) {
        const auto profile_iteration_start = std::chrono::steady_clock::now();
        const int obstacle_start = config.hpwl_only_iterations +
            config.progressive_density_iterations;
        const int stage = !config.progressive_legalization
            ? 0
            : iteration < config.hpwl_only_iterations ? 1
            : iteration < obstacle_start ? 2 : 3;
        if (stage != previous_stage) {
            if (previous_stage != 0) {
                std::fill(first.begin(), first.end(), 0.0);
                std::fill(second.begin(), second.end(), 0.0);
                std::fill(maximum_second.begin(), maximum_second.end(), 0.0);
                clear_bundle_models();
                integral_error = 0.0;
                serious_trial_scale = 1.0;
            }
            std::cout << "[Progressive] stage=" << stage
                      << " iter=" << iteration << '\n';
            previous_stage = stage;
        }
        dgx.clear(); dgy.clear(); fgx.clear(); fgy.clear();
        const auto profile_density_start = std::chrono::steady_clock::now();
        DensityResult d = density->compute(db, fillers, &dgx, &dgy, &fgx, &fgy);
        profile_density += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - profile_density_start).count();
        std::vector<Real> obstacle_gx, obstacle_gy;
        const ObstacleResult obstacle = obstacle_field
            ? obstacle_field->compute(
                  db, stage == 3 ? &obstacle_gx : nullptr,
                  stage == 3 ? &obstacle_gy : nullptr)
            : ObstacleResult{};
        if (config.multilevel_density &&
            (density->bins_x() < config.bins_x || density->bins_y() < config.bins_y)) {
            const bool coarse_stage = density->bins_x() < config.bins_x / 2 ||
                                      density->bins_y() < config.bins_y / 2;
            const int scheduled_iteration =
                (completed_grid_transitions + 1) * config.iterations /
                std::max(1, total_grid_transitions + 1);
            const bool advance = iteration >= scheduled_iteration ||
                                 (coarse_stage &&
                                  d.overflow <= config.multilevel_middle_overflow) ||
                                 (!coarse_stage &&
                                  d.overflow <= config.multilevel_fine_overflow);
            if (advance) {
                pending_effective_lambda = lambda_base * lambda_control;
                grid_x = std::min(config.bins_x, density->bins_x() * 2);
                grid_y = std::min(config.bins_y, density->bins_y() * 2);
                density = std::make_unique<ElectricDensity>(
                    db, grid_x, grid_y, config.target_density, field_model);
                d = density->compute(db, fillers, &dgx, &dgy, &fgx, &fgy);
                lambda_base = 0.0;
                ++completed_grid_transitions;
                last_grid_transition = iteration;
                clear_bundle_models();
                std::fill(first.begin(), first.end(), 0.0);
                std::fill(second.begin(), second.end(), 0.0);
                std::fill(maximum_second.begin(), maximum_second.end(), 0.0);
                serious_trial_scale = 1.0;
                // Overflow values from different grids are not comparable.
                // Rank checkpoints again on the new, finer discretization.
                result.have_feasible = false;
                best_feasible.clear();
                best_overflow = capture(db, fillers);
                best_overflow_metrics = Metrics{};
                best_overflow_metrics.overflow =
                    std::numeric_limits<Real>::infinity();
                best_legal_raw.clear();
                best_legal_hpwl = std::numeric_limits<Real>::infinity();
                std::cout << "[Multilevel] switched grid=" << grid_x << 'x'
                          << grid_y << " overflow=" << d.overflow << '\n';
            }
        }
        wgx.clear(); wgy.clear();
        std::vector<Real> group_hpwl;
        std::vector<std::vector<Real>> group_gx, group_gy;
        Real hpwl = 0.0;
        Real current_active_radius = 0.0;
        const bool grouped_bundle = config.enable_bundle && config.bundle_groups > 1;
        const auto profile_wire_start = std::chrono::steady_clock::now();
        if (grouped_bundle) {
            hpwl = exact_hpwl_group_subgradients(
                db, config.degree_limit, config.bundle_groups,
                group_hpwl, group_gx, group_gy);
            wgx.assign(db.nodes.size(), 0.0);
            wgy.assign(db.nodes.size(), 0.0);
            for (int group = 0; group < config.bundle_groups; ++group) {
                for (int id : db.movable_ids) {
                    wgx[id] += group_gx[group][id];
                    wgy[id] += group_gy[group][id];
                }
            }
        } else {
            Real active_radius = config.active_set_radius;
            if (refinement_active && config.refinement_active_set_radius >= 0.0) {
                const Real blend = config.refinement_active_set_decay_iterations > 0
                    ? std::clamp(
                          (iteration - refinement_start) / static_cast<Real>(
                              config.refinement_active_set_decay_iterations),
                          0.0, 1.0)
                    : 1.0;
                active_radius = (1.0 - blend) * config.active_set_radius +
                    blend * config.refinement_active_set_radius;
            }
            if (config.adaptive_active_set && active_radius > 0.0) {
                active_radius *= active_radius_scale;
            }
            current_active_radius = active_radius;
            hpwl = config.primal_dual_step > 0.0
                ? exact_hpwl_primal_dual_direction(
                      db, config.degree_limit, config.primal_dual_step,
                      &wgx, &wgy)
                : active_radius > 0.0
                ? exact_hpwl_active_set_direction(
                      db, config.degree_limit, active_radius,
                      config.active_set_power,
                      &wgx, &wgy,
                      config.adaptive_active_predictive
                          ? config.adaptive_active_set_span_cap : 0.0)
                : exact_hpwl_subgradient(
                      db, config.degree_limit, &wgx, &wgy);
        }
        profile_wire += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - profile_wire_start).count();
        ++result.objective_evaluations;
        const int active_sampling_count = refinement_active
            ? config.refinement_gradient_sampling_samples
            : config.gradient_sampling_samples;
        const Real active_sampling_radius = refinement_active
            ? config.refinement_gradient_sampling_radius
            : config.gradient_sampling_radius;
        if (active_sampling_count > 0 &&
            iteration % config.gradient_sampling_interval == 0) {
            sample_exact_subgradient(
                db, active_sampling_count,
                active_sampling_radius,
                config.seed ^ static_cast<std::uint64_t>(iteration + 1),
                config.degree_limit, wgx, wgy, result.objective_evaluations);
        }
        const bool density_active = iteration >= config.hpwl_only_iterations;
        if (density_active && lambda_base == 0.0) {
            Real wire_l1 = 0.0, density_l1 = 0.0;
            for (int id : db.movable_ids) {
                wire_l1 += std::abs(wgx[id]) + std::abs(wgy[id]);
                density_l1 += std::abs(dgx[id]) + std::abs(dgy[id]);
            }
            density_l1 += l1_norm(fgx) + l1_norm(fgy);
            lambda_base = config.density_weight_scale * wire_l1 /
                          std::max<Real>(density_l1, 1.0e-30);
            if (pending_effective_lambda > 0.0) {
                lambda_control = std::clamp(
                    pending_effective_lambda / std::max<Real>(lambda_base, 1.0e-30),
                    config.lambda_control_min, config.lambda_control_max);
                pending_effective_lambda = 0.0;
            }
            initial_density_overflow = d.overflow;
            previous_overflow = d.overflow;
            std::cout << "[Exact] lambda_base=" << lambda_base
                      << " wire_l1=" << wire_l1
                      << " density_l1=" << density_l1 << '\n';
        }
        const bool start_progressive_refinement =
            config.progressive_legalization && stage == 3 &&
            d.overflow <= config.progressive_overflow_upper;
        const Real refinement_start_overflow =
            config.refinement_start_overflow >= 0.0
                ? config.refinement_start_overflow : config.stop_overflow;
        const bool start_legacy_refinement = config.feasible_refinement &&
            !config.progressive_legalization &&
            d.overflow <= refinement_start_overflow;
        if (density_active && !refinement_active &&
            (start_progressive_refinement || start_legacy_refinement)) {
            refinement_active = true;
            refinement_start = iteration;
            if (config.adaptive_active_smart) {
                active_radius_scale = 1.0;
                adaptive_signal_ema = 0.0;
                adaptive_hpwl_history.clear();
                adaptive_overflow_history.clear();
            }
            std::fill(first.begin(), first.end(), 0.0);
            std::fill(second.begin(), second.end(), 0.0);
            std::fill(maximum_second.begin(), maximum_second.end(), 0.0);
            clear_bundle_models();
            integral_error = 0.0;
            serious_trial_scale = 1.0;
            std::cout << "[Refine] start iter=" << iteration
                      << " hpwl=" << hpwl << " overflow=" << d.overflow
                      << " optimizer=" << global_optimizer_name(config.refinement_optimizer)
                      << " lr_scale=" << config.refinement_learning_rate_scale << '\n';
        }
        const Real lambda = density_active ? lambda_base * lambda_control : 0.0;
        if (density_active && d.overflow <= config.bundle_start_overflow)
            add_row_proximity_force(db, config.legal_row_force, wgy);
        const Real progress = iteration / static_cast<Real>(std::max(1, config.iterations));
        const Real fraction = density_active ? config.nonsmooth_step_fraction
                                             : config.nonsmooth_hpwl_step_fraction;
        Real learning_rate = fraction * ((db.xh - db.xl) + (db.yh - db.yl)) /
                             std::sqrt(1.0 + 4.0 * progress);
        if (refinement_active) learning_rate *= config.refinement_learning_rate_scale;
        if (refinement_active && config.tangent_refinement)
            learning_rate *= tangent_trial_scale;
        const int transition_age = iteration - last_grid_transition;
        if (transition_age >= 0 && transition_age < 25) {
            learning_rate *= 0.10 + 0.90 * (transition_age + 1) / 25.0;
        }
        if (config.serious_bundle) learning_rate *= serious_trial_scale;
        const GlobalOptimizer active_optimizer = refinement_active
            ? config.refinement_optimizer : config.optimizer;
        BundleStats bundle_stats;
        if (config.enable_bundle && density_active && d.overflow <= config.bundle_start_overflow) {
            if (grouped_bundle) {
                const std::vector<Real> sampled_gx = wgx;
                const std::vector<Real> sampled_gy = wgy;
                std::vector<Real> grouped_gx, grouped_gy;
                bundle_stats = apply_grouped_bundle(
                    db, group_hpwl, learning_rate, iteration, config,
                    grouped_cuts, grouped_last_cut, bundle_prox_multiplier,
                    group_gx, group_gy, grouped_gx, grouped_gy);
                if (config.gradient_sampling_samples > 0) {
                    for (int id : db.movable_ids) {
                        wgx[id] = 0.5 * (sampled_gx[id] + grouped_gx[id]);
                        wgy[id] = 0.5 * (sampled_gy[id] + grouped_gy[id]);
                    }
                } else {
                    wgx = std::move(grouped_gx);
                    wgy = std::move(grouped_gy);
                }
            } else {
                bundle_stats = apply_bundle(db, hpwl, learning_rate, iteration, config,
                                            cuts, last_cut, bundle_prox_multiplier,
                                            wgx, wgy);
            }
        }

        Real obstacle_ramp = 0.0;
        Real obstacle_weight = 0.0;
        if (stage == 3) {
            const Real raw_ramp = (iteration - obstacle_start + 1) /
                static_cast<Real>(config.progressive_obstacle_iterations);
            obstacle_ramp = std::clamp(raw_ramp, 0.0, 1.0);
            obstacle_ramp = obstacle_ramp * obstacle_ramp *
                (3.0 - 2.0 * obstacle_ramp);
            if (obstacle_weight_base == 0.0) {
                Real wire_l1 = 0.0;
                Real obstacle_l1 = 0.0;
                for (int id : db.movable_ids) {
                    wire_l1 += std::abs(wgx[id]) + std::abs(wgy[id]);
                    obstacle_l1 += std::abs(obstacle_gx[id]) +
                                   std::abs(obstacle_gy[id]);
                }
                obstacle_weight_base = config.progressive_obstacle_scale *
                    wire_l1 / std::max<Real>(obstacle_l1, 1.0e-30);
                std::cout << "[Progressive] obstacle_base="
                          << obstacle_weight_base << " wire_l1=" << wire_l1
                          << " obstacle_l1=" << obstacle_l1 << '\n';
            }
            obstacle_weight = obstacle_ramp * obstacle_weight_base;
            std::vector<Real> legal_gx, legal_gy;
            compute_legalization_force(
                db, config.progressive_congestion_gain, legal_gx, legal_gy);
            const Real segment_strength = obstacle_ramp *
                config.progressive_segment_force;
            for (int id : db.movable_ids) {
                wgx[id] += segment_strength * legal_gx[id];
                wgy[id] += segment_strength * legal_gy[id];
            }
            add_row_proximity_force(
                db, obstacle_ramp * config.progressive_row_force, wgy);
        }

        LegalizationProxy legal_proxy;
        if (config.legal_checkpoint_selection || config.late_legal_projection ||
            config.legal_row_force > 0.0 || config.progressive_legalization) {
            legal_proxy = evaluate_legalization_proxy(db);
        }
        Real quick_legal_hpwl = 0.0;
        bool projection_applied = false;
        const bool checkpoint_due = config.legal_checkpoint_selection &&
            d.overflow <= config.stop_overflow &&
            iteration % config.legal_checkpoint_interval == 0;
        const bool projection_due = config.late_legal_projection &&
            d.overflow <= config.bundle_start_overflow &&
            iteration % config.legal_projection_interval == 0;
        if (checkpoint_due || projection_due) {
            Database legalized;
            if (quick_legalize(db, config.legal_checkpoint_detailed_passes,
                               legalized, quick_legal_hpwl)) {
                std::cout << "[QuickLegal] iter=" << iteration
                          << " gp_hpwl=" << hpwl
                          << " legal_hpwl=" << quick_legal_hpwl
                          << " row_distance=" << legal_proxy.mean_row_distance
                          << " segment_overflow=" << legal_proxy.segment_overflow << '\n';
                if (checkpoint_due && quick_legal_hpwl < best_legal_hpwl) {
                    best_legal_hpwl = quick_legal_hpwl;
                    best_legal_raw = capture(db, fillers);
                    best_legal_metrics.iteration = iteration;
                    best_legal_metrics.exact_hpwl = hpwl;
                    best_legal_metrics.overflow = d.overflow;
                }
                if (projection_due &&
                    quick_legal_hpwl <= hpwl * config.legal_projection_max_hpwl_ratio) {
                    const std::vector<Real> before_projection = capture(db, fillers);
                    const Real mix = config.legal_projection_mix;
                    for (int id : db.movable_ids) {
                        db.nodes[id].x = (1.0 - mix) * db.nodes[id].x +
                                         mix * legalized.nodes[id].x;
                        db.nodes[id].y = (1.0 - mix) * db.nodes[id].y +
                                         mix * legalized.nodes[id].y;
                    }
                    clamp_to_region(db, fillers);
                    const DensityResult projected_density = density->compute(
                        db, fillers, nullptr, nullptr, nullptr, nullptr);
                    const bool projection_density_safe = d.overflow <= config.stop_overflow
                        ? projected_density.overflow <= config.stop_overflow
                        : projected_density.overflow <= d.overflow + 0.002;
                    if (projection_density_safe) {
                        std::fill(first.begin(), first.end(), 0.0);
                        std::fill(second.begin(), second.end(), 0.0);
                        std::fill(maximum_second.begin(), maximum_second.end(), 0.0);
                        clear_bundle_models();
                        projection_applied = true;
                        std::cout << "[LegalProjection] iter=" << iteration
                                  << " mix=" << mix
                                  << " overflow=" << projected_density.overflow << '\n';
                    } else {
                        apply(db, fillers, before_projection);
                        std::cout << "[LegalProjection] iter=" << iteration
                                  << " rejected overflow="
                                  << projected_density.overflow << '\n';
                    }
                }
            }
        }

        Metrics current;
        current.iteration = iteration;
        current.exact_hpwl = hpwl;
        current.smooth_wirelength = hpwl;
        current.density_energy = d.energy;
        current.overflow = d.overflow;
        current.max_density = d.max_density;
        current.density_weight = lambda;
        current.step = learning_rate;
        if (d.overflow < best_overflow_metrics.overflow) {
            best_overflow_metrics = current;
            best_overflow = capture(db, fillers);
        }
        const Real feasible_overflow = config.progressive_legalization
            ? config.progressive_overflow_upper : config.stop_overflow;
        if (d.overflow <= feasible_overflow &&
            (!result.have_feasible || hpwl < result.best_feasible_metrics.exact_hpwl)) {
            result.have_feasible = true;
            result.best_feasible_metrics = current;
            best_feasible = capture(db, fillers);
        }
        if (stage == 3 && d.overflow <= config.progressive_overflow_upper) {
            const Real score = hpwl / config.progressive_hpwl_target +
                50.0 * obstacle.overlap_ratio +
                0.05 * legal_proxy.mean_row_distance +
                5.0 * legal_proxy.segment_overflow;
            if (score < best_progressive_score) {
                best_progressive_score = score;
                best_progressive = capture(db, fillers);
                best_progressive_metrics = current;
            }
        }
        if (metrics) {
            metrics << iteration << ',' << hpwl << ',' << d.overflow << ','
                    << d.max_density << ',' << d.energy << ',' << lambda_base << ','
                    << lambda_control << ',' << lambda << ',' << learning_rate << ','
                    << bundle_stats.cuts << ',' << bundle_stats.newest_weight << ','
                    << bundle_stats.model_error << ',' << bundle_stats.cosine << ','
                    << (refinement_active ? 1 : 0) << ','
                    << global_optimizer_name(active_optimizer) << ','
                    << density->bins_x() << ',' << density->bins_y() << ','
                    << legal_proxy.mean_row_distance << ','
                    << legal_proxy.segment_overflow << ','
                    << quick_legal_hpwl << ',' << last_serious_step << ','
                    << active_radius_scale << ',' << current_active_radius << ','
                    << adaptive_signal_ema << ',' << adaptive_phase;
            if (config.progressive_legalization) {
                metrics << ',' << stage << ',' << obstacle.energy << ','
                        << obstacle_weight << ',' << obstacle.overlap_area << ','
                        << obstacle.overlap_ratio << ',' << obstacle.overlapping_cells
                        << ',' << last_filter_fraction;
            }
            metrics << '\n';
        }
        if (iteration % config.log_every == 0 || iteration + 1 == config.iterations) {
            std::cout << "[Exact] iter=" << iteration << " hpwl=" << hpwl
                      << " overflow=" << d.overflow << " eta=" << lambda_control
                      << " lambda=" << lambda << " lr=" << learning_rate
                      << " cuts=" << bundle_stats.cuts
                      << " grid=" << density->bins_x() << 'x' << density->bins_y()
                      << " refine=" << (refinement_active ? 1 : 0)
                      << " stage=" << stage
                      << " macro_overlap=" << obstacle.overlap_ratio
                      << " macro_cells=" << obstacle.overlapping_cells << '\n';
        }
        save_global_snapshot(db, config, iteration);

        if (projection_applied) {
            previous_hpwl = hpwl;
            continue;
        }

        std::fill(gradient.begin(), gradient.end(), 0.0);
        if (refinement_active && config.tangent_refinement) {
            std::fill(wire_gradient_refinement.begin(),
                      wire_gradient_refinement.end(), 0.0);
            std::fill(density_normal_refinement.begin(),
                      density_normal_refinement.end(), 0.0);
            std::vector<Real>& wire_gradient = wire_gradient_refinement;
            std::vector<Real>& density_normal = density_normal_refinement;
            const int movable_count = static_cast<int>(db.movable_ids.size());
            const int filler_count = static_cast<int>(fillers.size());
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < movable_count; ++i) {
                const int id = db.movable_ids[i];
                const Node& node = db.nodes[id];
                wire_gradient[i] = wgx[id] /
                    std::max<Real>(1.0, db.node_pin_weight[id]);
                density_normal[i] = dgx[id] /
                    std::max<Real>(1.0, node.area());
                wire_gradient[n + i] = wgy[id] /
                    std::max<Real>(1.0, db.node_pin_weight[id]);
                density_normal[n + i] = dgy[id] /
                    std::max<Real>(1.0, node.area());
            }
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < filler_count; ++i) {
                const std::size_t index = db.movable_ids.size() + i;
                const Real area = std::max<Real>(1.0, fillers[i].width * fillers[i].height);
                density_normal[index] = fgx[i] / area;
                density_normal[n + index] = fgy[i] / area;
            }

            Real wire_density_dot = 0.0;
            Real density_norm2 = 0.0;
            Real wire_norm2 = 0.0;
            #pragma omp parallel for reduction(+:wire_density_dot,density_norm2,wire_norm2) schedule(static)
            for (int i = 0; i < static_cast<int>(gradient.size()); ++i) {
                wire_density_dot += wire_gradient[i] * density_normal[i];
                density_norm2 += density_normal[i] * density_normal[i];
                wire_norm2 += wire_gradient[i] * wire_gradient[i];
            }
            const Real tangent_coefficient = wire_density_dot < 0.0
                ? wire_density_dot / std::max<Real>(density_norm2, 1.0e-30)
                : 0.0;
            const Real overflow_band = std::max<Real>(
                config.stop_overflow - config.refinement_lower_overflow, 1.0e-4);
            const Real recovery_fraction = std::clamp(
                (d.overflow - config.stop_overflow) / overflow_band, 0.0, 1.0);
            const Real recovery_scale = recovery_fraction * std::sqrt(
                wire_norm2 / std::max<Real>(density_norm2, 1.0e-30));
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < static_cast<int>(gradient.size()); ++i) {
                gradient[i] = wire_gradient[i] -
                    tangent_coefficient * density_normal[i] +
                    recovery_scale * density_normal[i];
            }
        } else {
            const int movable_count = static_cast<int>(db.movable_ids.size());
            const int filler_count = static_cast<int>(fillers.size());
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < movable_count; ++i) {
                const int id = db.movable_ids[i];
                const Node& node = db.nodes[id];
                const Real preconditioner = std::max<Real>(
                    1.0, db.node_pin_weight[id] + lambda * node.area());
                gradient[i] = (wgx[id] + lambda * dgx[id] +
                                 obstacle_weight * (stage == 3 ? obstacle_gx[id] : 0.0)) /
                                preconditioner;
                gradient[n + i] = (wgy[id] + lambda * dgy[id] +
                                   obstacle_weight * (stage == 3 ? obstacle_gy[id] : 0.0)) /
                                  preconditioner;
            }
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < filler_count; ++i) {
                const std::size_t index = db.movable_ids.size() + i;
                const Real preconditioner = std::max<Real>(
                    1.0, lambda * fillers[i].width * fillers[i].height);
                gradient[index] = lambda * fgx[i] / preconditioner;
                gradient[n + index] = lambda * fgy[i] / preconditioner;
            }
        }

        capture_into(db, fillers, previous_positions);
        positions = previous_positions;
        const int age = refinement_active ? iteration - refinement_start + 1
                                          : iteration + 1;
        const Real first_bias_denominator = std::max<Real>(
            1.0e-12, 1.0 - std::pow(config.beta1, age));
        const Real second_bias_denominator = std::max<Real>(
            1.0e-12, 1.0 - std::pow(config.beta2, age));
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(positions.size()); ++i) {
            positions[i] -= coordinate_step(active_optimizer, gradient[i], learning_rate,
                                             config, first[i], second[i],
                                             maximum_second[i],
                                             first_bias_denominator,
                                             second_bias_denominator);
        }
        apply(db, fillers, positions);
        if (config.progressive_legalization && config.progressive_filter) {
            const std::vector<Real> full_trial = positions;
            bool accepted = false;
            Real fraction = 1.0;
            for (int backtrack = 0;
                 backtrack <= config.progressive_filter_backtracks;
                 ++backtrack, fraction *= 0.5) {
                std::vector<Real> trial(previous_positions.size());
                for (std::size_t i = 0; i < trial.size(); ++i) {
                    trial[i] = previous_positions[i] +
                        fraction * (full_trial[i] - previous_positions[i]);
                }
                apply(db, fillers, trial);
                const Real candidate_hpwl = exact_hpwl(db);
                ++result.objective_evaluations;
                if (stage == 1) {
                    if (candidate_hpwl <= hpwl + 1.0e-9) {
                        positions = std::move(trial);
                        accepted = true;
                        break;
                    }
                    continue;
                }
                if (stage == 2 &&
                    d.overflow > config.progressive_overflow_upper + 0.02) {
                    positions = std::move(trial);
                    accepted = true;
                    break;
                }
                const DensityResult candidate_density = density->compute(
                    db, fillers, nullptr, nullptr, nullptr, nullptr);
                const ObstacleResult candidate_obstacle = obstacle_field->compute(
                    db, nullptr, nullptr);
                const Real target = config.progressive_hpwl_target;
                const Real upper = config.progressive_overflow_upper;
                const Real current_overflow_violation = std::max<Real>(
                    0.0, (d.overflow - upper) / std::max<Real>(upper, 1.0e-6));
                const Real candidate_overflow_violation = std::max<Real>(
                    0.0, (candidate_density.overflow - upper) /
                         std::max<Real>(upper, 1.0e-6));
                const Real current_merit = hpwl / target +
                    4.0 * current_overflow_violation +
                    100.0 * obstacle_ramp * obstacle.overlap_ratio;
                const Real candidate_merit = candidate_hpwl / target +
                    4.0 * candidate_overflow_violation +
                    100.0 * obstacle_ramp * candidate_obstacle.overlap_ratio;
                const Real dynamic_hpwl_ratio = 1.067 + 0.40 * std::min<Real>(
                    1.0, current_overflow_violation);
                const bool hpwl_safe = d.overflow > upper + 0.002 ||
                                       candidate_hpwl <=
                                           dynamic_hpwl_ratio * target ||
                                       candidate_hpwl <= hpwl;
                const bool overflow_safe = d.overflow <= upper + 0.002
                    ? candidate_density.overflow <= upper + 0.003
                    : candidate_density.overflow <= d.overflow + 0.005;
                const Real obstacle_slack = stage == 3
                    ? std::max<Real>(1.0e-8, (1.0 - obstacle_ramp) * 1.0e-4)
                    : std::numeric_limits<Real>::infinity();
                const bool obstacle_safe = candidate_obstacle.overlap_ratio <=
                    obstacle.overlap_ratio + obstacle_slack;
                if (hpwl_safe && overflow_safe && obstacle_safe &&
                    candidate_merit <= current_merit + 1.0e-4) {
                    positions = std::move(trial);
                    accepted = true;
                    break;
                }
            }
            if (!accepted) {
                positions = previous_positions;
                apply(db, fillers, positions);
                std::fill(first.begin(), first.end(), 0.0);
                std::fill(second.begin(), second.end(), 0.0);
                std::fill(maximum_second.begin(), maximum_second.end(), 0.0);
                last_filter_fraction = 0.0;
            } else {
                apply(db, fillers, positions);
                last_filter_fraction = fraction;
            }
        }
        if (config.tangent_refinement && refinement_active &&
            !config.progressive_legalization) {
            const std::vector<Real> full_trial = positions;
            bool accepted = false;
            Real accepted_fraction = 0.0;
            for (int backtrack = 0;
                 backtrack <= config.refinement_filter_backtracks;
                 ++backtrack) {
                const Real fraction = std::ldexp(1.0, -backtrack);
                std::vector<Real> trial(previous_positions.size());
                #pragma omp parallel for schedule(static)
                for (int i = 0; i < static_cast<int>(trial.size()); ++i) {
                    trial[i] = previous_positions[i] +
                        fraction * (full_trial[i] - previous_positions[i]);
                }
                apply(db, fillers, trial);
                const Real candidate_hpwl = exact_hpwl(db);
                ++result.objective_evaluations;
                const DensityResult candidate_density = density->compute(
                    db, fillers, nullptr, nullptr, nullptr, nullptr);
                if (candidate_hpwl <= hpwl + 1.0e-9 &&
                    candidate_density.overflow <= config.stop_overflow) {
                    positions = std::move(trial);
                    accepted = true;
                    accepted_fraction = fraction;
                    break;
                }
            }
            if (accepted) {
                apply(db, fillers, positions);
                tangent_trial_scale = std::min<Real>(
                    1.0, 1.05 * tangent_trial_scale);
            } else {
                apply(db, fillers, previous_positions);
                std::fill(first.begin(), first.end(), 0.0);
                std::fill(second.begin(), second.end(), 0.0);
                std::fill(maximum_second.begin(), maximum_second.end(), 0.0);
                tangent_trial_scale = std::max<Real>(
                    1.0e-3, 0.5 * tangent_trial_scale);
            }
            if ((!accepted || accepted_fraction < 1.0) &&
                iteration % config.log_every == 0) {
                std::cout << "[TangentFilter] iter=" << iteration
                          << " accepted=" << (accepted ? 1 : 0)
                          << " fraction=" << accepted_fraction
                          << " trial_scale=" << tangent_trial_scale << '\n';
            }
        }
        last_serious_step = 1;
        if (config.serious_bundle && bundle_stats.cuts > 0 &&
            d.overflow <= config.bundle_start_overflow) {
            std::vector<Real> trial_group_hpwl;
            std::vector<std::vector<Real>> trial_group_gx, trial_group_gy;
            Real candidate_hpwl = grouped_bundle
                ? exact_hpwl_group_subgradients(
                      db, config.degree_limit, config.bundle_groups,
                      trial_group_hpwl, trial_group_gx, trial_group_gy)
                : exact_hpwl(db);
            ++result.objective_evaluations;
            DensityResult candidate_density = density->compute(
                db, fillers, nullptr, nullptr, nullptr, nullptr);
            ObstacleResult candidate_obstacle = obstacle_field
                ? obstacle_field->compute(db, nullptr, nullptr)
                : ObstacleResult{};
            const Real bundle_feasible_overflow = config.progressive_legalization
                ? config.progressive_overflow_upper : config.stop_overflow;
            if (d.overflow <= bundle_feasible_overflow &&
                candidate_density.overflow > bundle_feasible_overflow &&
                candidate_density.overflow <=
                    bundle_feasible_overflow + config.bundle_overflow_tolerance) {
                const std::vector<Real> unprojected = positions;
                std::vector<Real> feasible = previous_positions;
                DensityResult feasible_density = d;
                Real lower = 0.0, upper = 1.0;
                for (int search = 0; search < 10; ++search) {
                    const Real fraction = 0.5 * (lower + upper);
                    std::vector<Real> trial(previous_positions.size());
                    for (std::size_t i = 0; i < trial.size(); ++i) {
                        trial[i] = previous_positions[i] +
                            fraction * (unprojected[i] - previous_positions[i]);
                    }
                    apply(db, fillers, trial);
                    const DensityResult trial_density = density->compute(
                        db, fillers, nullptr, nullptr, nullptr, nullptr);
                    if (trial_density.overflow <= bundle_feasible_overflow) {
                        lower = fraction;
                        feasible = std::move(trial);
                        feasible_density = trial_density;
                    } else {
                        upper = fraction;
                    }
                }
                positions = std::move(feasible);
                apply(db, fillers, positions);
                candidate_density = feasible_density;
                if (obstacle_field) {
                    candidate_obstacle = obstacle_field->compute(
                        db, nullptr, nullptr);
                }
                if (grouped_bundle) {
                    candidate_hpwl = exact_hpwl_group_subgradients(
                        db, config.degree_limit, config.bundle_groups,
                        trial_group_hpwl, trial_group_gx, trial_group_gy);
                } else {
                    candidate_hpwl = exact_hpwl(db);
                }
                ++result.objective_evaluations;
                std::cout << "[BundleBoundary] iter=" << iteration
                          << " fraction=" << lower
                          << " hpwl=" << candidate_hpwl
                          << " overflow=" << candidate_density.overflow << '\n';
            }
            Real predicted_decrease = 0.0;
            for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
                const int id = db.movable_ids[i];
                predicted_decrease -= (wgx[id] + lambda * dgx[id]) *
                    (positions[i] - previous_positions[i]);
                predicted_decrease -= (wgy[id] + lambda * dgy[id]) *
                    (positions[n + i] - previous_positions[n + i]);
                if (stage == 3) {
                    predicted_decrease -= obstacle_weight * obstacle_gx[id] *
                        (positions[i] - previous_positions[i]);
                    predicted_decrease -= obstacle_weight * obstacle_gy[id] *
                        (positions[n + i] - previous_positions[n + i]);
                }
            }
            for (std::size_t i = 0; i < fillers.size(); ++i) {
                const std::size_t index = db.movable_ids.size() + i;
                predicted_decrease -= lambda * fgx[i] *
                    (positions[index] - previous_positions[index]);
                predicted_decrease -= lambda * fgy[i] *
                    (positions[n + index] - previous_positions[n + index]);
            }
            const bool in_filter_band = d.overflow <=
                bundle_feasible_overflow + config.bundle_overflow_tolerance;
            const Real actual_decrease = config.progressive_legalization
                ? (hpwl + lambda * d.energy + obstacle_weight * obstacle.energy) -
                  (candidate_hpwl + lambda * candidate_density.energy +
                   obstacle_weight * candidate_obstacle.energy)
                : in_filter_band
                    ? hpwl - candidate_hpwl
                    : (hpwl + lambda * d.energy) -
                      (candidate_hpwl + lambda * candidate_density.energy);
            // Grid overflow is nonconvex and can rise slightly even when the
            // electrostatic merit falls.  Permit a small infeasible-stage
            // trust-band excursion, but never leave the feasible region once
            // the requested overflow threshold has been reached.
            constexpr Real infeasible_overflow_tolerance = 0.002;
            const bool density_safe = in_filter_band
                ? candidate_density.overflow <=
                      bundle_feasible_overflow + config.bundle_overflow_tolerance
                : candidate_density.overflow <=
                      d.overflow + infeasible_overflow_tolerance;
            const bool serious = density_safe && actual_decrease > 0.0 &&
                (predicted_decrease <= 1.0e-12 ||
                 actual_decrease >= config.serious_step_ratio * predicted_decrease);
            if (serious) {
                bundle_prox_multiplier = std::min<Real>(4.0,
                                                        1.10 * bundle_prox_multiplier);
                serious_trial_scale = std::min<Real>(1.0,
                                                      1.15 * serious_trial_scale);
                std::cout << "[BundleSerious] iter=" << iteration
                          << " actual=" << actual_decrease
                          << " predicted=" << predicted_decrease
                          << " prox=" << bundle_prox_multiplier
                          << " trial_scale=" << serious_trial_scale << '\n';
            } else {
                if (grouped_bundle) {
                    GlobalPlaceConfig group_config = config;
                    group_config.bundle_size = config.bundle_group_size;
                    for (int group = 0; group < config.bundle_groups; ++group) {
                        append_bundle_cut(
                            db, trial_group_hpwl[group], trial_group_gx[group],
                            trial_group_gy[group], group_config,
                            grouped_cuts[group]);
                    }
                } else {
                    std::vector<Real> trial_gx, trial_gy;
                    exact_hpwl_subgradient(db, config.degree_limit,
                                           &trial_gx, &trial_gy);
                    ++result.objective_evaluations;
                    append_bundle_cut(db, candidate_hpwl, trial_gx, trial_gy,
                                      config, cuts);
                }
                apply(db, fillers, previous_positions);
                std::fill(first.begin(), first.end(), 0.0);
                std::fill(second.begin(), second.end(), 0.0);
                std::fill(maximum_second.begin(), maximum_second.end(), 0.0);
                bundle_prox_multiplier = std::max<Real>(0.10,
                                                        0.50 * bundle_prox_multiplier);
                serious_trial_scale = std::max<Real>(0.01,
                                                     0.50 * serious_trial_scale);
                last_serious_step = 0;
                std::cout << "[BundleNull] iter=" << iteration
                          << " actual=" << actual_decrease
                          << " predicted=" << predicted_decrease
                          << " overflow_trial=" << candidate_density.overflow
                          << " density_safe=" << (density_safe ? 1 : 0)
                          << " prox=" << bundle_prox_multiplier
                          << " trial_scale=" << serious_trial_scale << '\n';
            }
        }
        if (density_active) {
            ++density_step;
            if (config.progressive_legalization) {
                if (d.overflow > config.progressive_dual_switch_overflow) {
                    lambda_control = update_lambda_control(
                        lambda_control, LambdaPolicy::Dreamplace, density_step,
                        hpwl, d.overflow, previous_hpwl,
                        initial_density_overflow, previous_overflow,
                        integral_error, config);
                } else {
                    lambda_control = progressive_band_dual_lambda(
                        lambda_control, density_step, hpwl, d.overflow,
                        initial_density_overflow, integral_error, config);
                }
            } else if (refinement_active) {
                lambda_control = feasible_band_lambda(lambda_control, d.overflow, config);
            } else {
                lambda_control = update_lambda_control(
                    lambda_control, config.lambda_policy, density_step, hpwl, d.overflow,
                    previous_hpwl, initial_density_overflow, previous_overflow,
                    integral_error, config);
            }
        }
        adaptive_hpwl_history.push_back(hpwl);
        adaptive_overflow_history.push_back(d.overflow);
        while (adaptive_hpwl_history.size() > static_cast<std::size_t>(
                   config.adaptive_active_set_window)) {
            adaptive_hpwl_history.pop_front();
            adaptive_overflow_history.pop_front();
        }
        if (config.adaptive_active_set && config.active_set_radius > 0.0 &&
            (!refinement_active || config.adaptive_active_set_refinement) &&
            (iteration + 1) % config.adaptive_active_set_interval == 0) {
            const Real upper = config.stop_overflow + 0.010;
            const Real lower = std::max<Real>(0.0, config.refinement_lower_overflow);
            const Real hpwl_rise = hpwl > previous_hpwl * 1.001;
            const Real old_scale = active_radius_scale;
            if (config.adaptive_active_smart && adaptive_hpwl_history.size() >=
                static_cast<std::size_t>(std::max(4, config.adaptive_active_set_window / 2))) {
                // Use a dead-banded, exponentially filtered signal.  Overflow
                // is the hard constraint; HPWL only chooses the radius inside
                // the feasible band.  The log-step cap prevents oscillations
                // from changing the near-extreme set too abruptly.
                const std::size_t half = std::max<std::size_t>(
                    2, adaptive_hpwl_history.size() / 2);
                const std::size_t split = adaptive_hpwl_history.size() - half;
                const Real hpwl_old = adaptive_hpwl_history[split - 1];
                const Real overflow_old = adaptive_overflow_history[split - 1];
                const Real hpwl_trend = (hpwl - hpwl_old) /
                    std::max<Real>(std::abs(hpwl_old), 1.0);
                const Real overflow_trend = d.overflow - overflow_old;
                const Real band = std::max<Real>(
                    config.stop_overflow - config.refinement_lower_overflow, 0.005);
                const Real deadband = std::max<Real>(
                    config.adaptive_active_set_deadband, 0.0);
                const Real pressure = (d.overflow - config.stop_overflow) / band;
                const Real trend = overflow_trend / band;
                Real signal = 0.0;
                if (config.adaptive_active_predictive) {
                    const Real far_threshold = std::max<Real>(
                        config.stop_overflow + 0.08, 2.0 * config.stop_overflow);
                    if (refinement_active) adaptive_phase = 3;
                    else if (d.overflow > far_threshold) adaptive_phase = 1;
                    else adaptive_phase = 2;

                    const Real expected_drop = std::max<Real>(
                        0.0, (overflow_old - config.stop_overflow) *
                        static_cast<Real>(half) /
                        std::max<Real>(1.0, config.lambda_trajectory_horizon));
                    const Real actual_drop = overflow_old - d.overflow;
                    const Real lag = (expected_drop - actual_drop) / band;
                    const Real hpwl_term = std::clamp(
                        hpwl_trend / 0.01, -1.0, 1.0);
                    if (adaptive_phase == 1) {
                        signal = 0.55 * std::clamp(lag, -1.0, 2.0) +
                                 0.10 * std::clamp(pressure, 0.0, 2.0) -
                                 0.10 * hpwl_term;
                    } else if (adaptive_phase == 2) {
                        signal = 0.40 * std::clamp(pressure, -1.0, 1.5) +
                                 0.25 * std::clamp(lag, -1.0, 1.0) -
                                 0.45 * hpwl_term;
                    } else {
                        signal = 0.25 * std::clamp(pressure, -1.0, 1.0) -
                                 0.65 * hpwl_term;
                    }
                } else {
                    if (d.overflow > config.stop_overflow + deadband)
                        signal += 0.75 * pressure + 0.25 * std::max<Real>(trend, 0.0);
                    else if (d.overflow < config.refinement_lower_overflow - deadband)
                        signal -= 0.65 * (-pressure) + 0.20 * std::max<Real>(-trend, 0.0);
                    else if (hpwl_trend > 0.0015)
                        signal -= 0.55 * std::min<Real>(hpwl_trend / 0.01, 1.0);
                    else if (hpwl_trend < -0.0015)
                        signal += 0.20 * std::min<Real>(-hpwl_trend / 0.01, 1.0);
                }
                const Real ema_keep = config.adaptive_active_predictive ? 0.80 : 0.70;
                adaptive_signal_ema = ema_keep * adaptive_signal_ema +
                    (1.0 - ema_keep) * signal;
                const Real configured_max_step = refinement_active &&
                    config.adaptive_active_predictive
                    ? config.adaptive_active_set_refinement_max_log_step
                    : config.adaptive_active_set_max_log_step;
                const Real max_log_step = std::max<Real>(configured_max_step, 0.002);
                const Real log_step = std::clamp(
                    0.07 * config.adaptive_active_set_gain * adaptive_signal_ema,
                    -max_log_step, max_log_step);
                active_radius_scale *= std::exp(log_step);
                std::cout << "[AdaptiveEpsilonSmart] iter=" << iteration
                          << " hpwl_trend=" << hpwl_trend
                          << " overflow_trend=" << overflow_trend
                          << " signal=" << adaptive_signal_ema
                          << " phase=" << adaptive_phase
                          << " log_step=" << log_step;
            } else if (d.overflow > upper) {
                active_radius_scale *= 1.10;
            } else if (d.overflow < lower) {
                active_radius_scale *= 0.90;
            } else if (hpwl_rise && d.overflow <= config.stop_overflow) {
                active_radius_scale *= 0.85;
            } else if (!hpwl_rise && d.overflow <= config.stop_overflow + 0.002) {
                active_radius_scale *= 1.03;
            }
            active_radius_scale = std::clamp(
                active_radius_scale,
                config.adaptive_active_set_min_scale,
                config.adaptive_active_set_max_scale);
            if (std::abs(active_radius_scale - old_scale) > 1.0e-12) {
                std::cout << "[AdaptiveEpsilon] iter=" << iteration
                          << " overflow=" << d.overflow
                          << " hpwl=" << hpwl
                          << " scale=" << active_radius_scale << '\n';
            }
        }
        profile_update += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - profile_iteration_start).count();
        previous_hpwl = hpwl;
    }

    if (config.profile) {
        std::cout << "[Profile] density_seconds=" << profile_density
                  << " wirelength_seconds=" << profile_wire
                  << " total_loop_seconds=" << profile_update
                  << " other_seconds=" << std::max<Real>(
                      0.0, profile_update - profile_density - profile_wire)
                  << '\n';
    }

    if (config.progressive_legalization && !best_progressive.empty()) {
        apply(db, fillers, best_progressive);
        result.have_feasible = true;
        result.best_feasible_metrics = best_progressive_metrics;
        std::cout << "[Progressive] restored stage-3 score="
                  << best_progressive_score << " iter="
                  << best_progressive_metrics.iteration << " hpwl="
                  << best_progressive_metrics.exact_hpwl << " overflow="
                  << best_progressive_metrics.overflow << '\n';
    } else if (config.legal_checkpoint_selection && !best_legal_raw.empty()) {
        apply(db, fillers, best_legal_raw);
        result.have_feasible = true;
        result.best_feasible_metrics = best_legal_metrics;
        result.selected_legal_hpwl = best_legal_hpwl;
        result.selected_legal_iteration = best_legal_metrics.iteration;
        std::cout << "[Exact] restored legal-aware iter="
                  << best_legal_metrics.iteration
                  << " gp_hpwl=" << best_legal_metrics.exact_hpwl
                  << " overflow=" << best_legal_metrics.overflow
                  << " quick_legal_hpwl=" << best_legal_hpwl << '\n';
    } else if (result.have_feasible) {
        apply(db, fillers, best_feasible);
        std::cout << "[Exact] restored feasible iter="
                  << result.best_feasible_metrics.iteration << " hpwl="
                  << result.best_feasible_metrics.exact_hpwl << " overflow="
                  << result.best_feasible_metrics.overflow << '\n';
    } else {
        apply(db, fillers, best_overflow);
        std::cout << "[Exact] no feasible state; restored min-overflow iter="
                  << best_overflow_metrics.iteration << " overflow="
                  << best_overflow_metrics.overflow << '\n';
    }
    const DensityResult final_density = density->compute(
        db, fillers, nullptr, nullptr, nullptr, nullptr);
    result.final_metrics.exact_hpwl = exact_hpwl(db);
    result.final_metrics.smooth_wirelength = result.final_metrics.exact_hpwl;
    result.final_metrics.overflow = final_density.overflow;
    result.final_metrics.max_density = final_density.max_density;
    result.final_metrics.density_weight = lambda_base * lambda_control;
    result.wall_time_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return result;
}

}  // namespace dpcpp
