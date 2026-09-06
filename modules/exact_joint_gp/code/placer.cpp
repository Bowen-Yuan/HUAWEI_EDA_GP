#include "epsilon_active/placer.hpp"

#include "epsilon_active/bookshelf.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <vector>

#include <omp.h>

namespace ea {
namespace {

std::vector<Real> capture(const Database& db) {
    const std::size_t n = db.movable_ids.size();
    std::vector<Real> positions(2 * n);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const Node& node = db.nodes[db.movable_ids[i]];
        positions[i] = node.x;
        positions[n + i] = node.y;
    }
    return positions;
}

void apply_positions(Database& db, const std::vector<Real>& positions) {
    const std::size_t n = db.movable_ids.size();
    if (positions.size() != 2 * n) {
        throw std::invalid_argument("position vector has wrong size");
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        Node& node = db.nodes[db.movable_ids[i]];
        node.x = std::clamp(positions[i], db.xl + 0.5 * node.width,
                            db.xh - 0.5 * node.width);
        node.y = std::clamp(positions[n + i], db.yl + 0.5 * node.height,
                            db.yh - 0.5 * node.height);
    }
}

void save_snapshot(const Database& db, const PlaceConfig& config, int iteration) {
    if (config.snapshot_every <= 0 || iteration % config.snapshot_every != 0) return;
    std::ostringstream name;
    name << "global_" << std::setw(5) << std::setfill('0') << iteration << ".pl";
    write_bookshelf_placement(db, config.output_dir / "snapshots" / name.str());
}

void configure_threads(int requested) {
    if (requested < 0 || requested > 40) {
        throw std::invalid_argument("thread count must be in [0, 40]");
    }
    omp_set_dynamic(0);
    const int threads = requested == 0 ? std::min(40, omp_get_max_threads()) : requested;
    omp_set_num_threads(std::max(1, threads));
}

void compute_net_batch_direction(const Database& db,
                                 const std::vector<Real>& node_x,
                                 const std::vector<Real>& node_y,
                                 int degree_limit,
                                 std::vector<Real>& batch_x,
                                 std::vector<Real>& batch_y) {
    batch_x.assign(db.nodes.size(), 0.0);
    batch_y.assign(db.nodes.size(), 0.0);
    std::vector<int> counts(db.nodes.size(), 0);
    for (const Net& net : db.nets) {
        if (net.pin_count < 2 ||
            static_cast<int>(net.pin_count) > degree_limit) continue;
        Real mean_x = 0.0;
        Real mean_y = 0.0;
        int count = 0;
        for (std::size_t p = net.pin_begin;
             p < net.pin_begin + net.pin_count; ++p) {
            const int id = db.pins[p].node;
            if (db.nodes[id].fixed) continue;
            mean_x += node_x[id];
            mean_y += node_y[id];
            ++count;
        }
        if (count == 0) continue;
        mean_x /= static_cast<Real>(count);
        mean_y /= static_cast<Real>(count);
        for (std::size_t p = net.pin_begin;
             p < net.pin_begin + net.pin_count; ++p) {
            const int id = db.pins[p].node;
            if (db.nodes[id].fixed) continue;
            batch_x[id] += mean_x;
            batch_y[id] += mean_y;
            ++counts[id];
        }
    }
    #pragma omp parallel for schedule(static)
    for (int id = 0; id < static_cast<int>(db.nodes.size()); ++id) {
        if (counts[id] > 0) {
            batch_x[id] /= static_cast<Real>(counts[id]);
            batch_y[id] /= static_cast<Real>(counts[id]);
        }
    }
}

void share_density_direction_over_nets(
    const Database& db, const DensityNetShareConfig& config,
    std::vector<Real>& density_x, std::vector<Real>& density_y) {
    if (config.weight <= 0.0) return;
    const Real blend = std::clamp(config.weight, 0.0, 1.0);
    std::vector<Real> shared_x(db.nodes.size(), 0.0);
    std::vector<Real> shared_y(db.nodes.size(), 0.0);
    std::vector<int> counts(db.nodes.size(), 0);
    for (const Net& net : db.nets) {
        if (net.pin_count < 2 ||
            static_cast<int>(net.pin_count) > config.degree_limit) continue;
        Real mean_x = 0.0;
        Real mean_y = 0.0;
        int movable = 0;
        for (std::size_t p = net.pin_begin;
             p < net.pin_begin + net.pin_count; ++p) {
            const int id = db.pins[p].node;
            if (db.nodes[id].fixed) continue;
            mean_x += density_x[id];
            mean_y += density_y[id];
            ++movable;
        }
        if (movable == 0) continue;
        mean_x /= static_cast<Real>(movable);
        mean_y /= static_cast<Real>(movable);
        for (std::size_t p = net.pin_begin;
             p < net.pin_begin + net.pin_count; ++p) {
            const int id = db.pins[p].node;
            if (db.nodes[id].fixed) continue;
            shared_x[id] += mean_x;
            shared_y[id] += mean_y;
            ++counts[id];
        }
    }
    #pragma omp parallel for schedule(static)
    for (int id = 0; id < static_cast<int>(db.nodes.size()); ++id) {
        if (counts[id] == 0) continue;
        const Real sx = shared_x[id] / static_cast<Real>(counts[id]);
        const Real sy = shared_y[id] / static_cast<Real>(counts[id]);
        const Real norm = std::hypot(density_x[id], density_y[id]);
        // For zero-gradient nodes, sharing is the only available exact
        // active-set direction.  For nonzero nodes, use a conservative blend.
        const Real local_blend = norm <= config.zero_gradient_threshold
            ? 1.0 : blend;
        density_x[id] = (1.0 - local_blend) * density_x[id] + local_blend * sx;
        density_y[id] = (1.0 - local_blend) * density_y[id] + local_blend * sy;
    }
}

void project_density_direction_from_hpwl(
    const Database& db, Real projection,
    const std::vector<Real>& wire_x, const std::vector<Real>& wire_y,
    std::vector<Real>& density_x, std::vector<Real>& density_y) {
    if (projection <= 0.0) return;
    const Real alpha = std::clamp(projection, 0.0, 1.0);
    #pragma omp parallel for schedule(static)
    for (int id = 0; id < static_cast<int>(db.nodes.size()); ++id) {
        const Real wx = wire_x[id];
        const Real wy = wire_y[id];
        const Real norm2 = wx * wx + wy * wy;
        if (norm2 <= 1.0e-24) continue;
        const Real coefficient = alpha *
            (density_x[id] * wx + density_y[id] * wy) / norm2;
        density_x[id] -= coefficient * wx;
        density_y[id] -= coefficient * wy;
    }
}

void update_regional_prices(const ExactOverlapDensity& density,
                            const RegionalPriceConfig& config,
                            std::vector<Real>& prices) {
    if (!config.enabled) return;
    const auto& occupancy = density.occupancy();
    const Real capacity = density.target_density() * density.bin_area();
    if (capacity <= 0.0 || occupancy.size() != prices.size()) return;
    Real mean = 0.0;
    int active = 0;
    for (std::size_t i = 0; i < occupancy.size(); ++i) {
        const Real ratio = std::max<Real>(
            (occupancy[i] - capacity) / capacity, 0.0);
        prices[i] = std::clamp(
            prices[i] + config.rho * (ratio - config.target_excess),
            config.min_price, config.max_price);
        if (ratio > 0.0) {
            mean += prices[i];
            ++active;
        }
    }
    if (active > 0) {
        mean /= static_cast<Real>(active);
        if (mean > 0.0) {
            for (std::size_t i = 0; i < prices.size(); ++i) {
                prices[i] = std::clamp(prices[i] / mean,
                                       config.min_price, config.max_price);
            }
        }
    }
}

class AdaptiveEpsilon {
public:
    explicit AdaptiveEpsilon(const PlaceConfig& config) : config_(config) {}

    void observe(int iteration, Real hpwl, Real overflow) {
        hpwl_history_.push_back(hpwl);
        overflow_history_.push_back(overflow);
        while (hpwl_history_.size() > static_cast<std::size_t>(config_.adaptive_window)) {
            hpwl_history_.pop_front();
            overflow_history_.pop_front();
        }
        if (!config_.adaptive_epsilon || (iteration + 1) % config_.adaptive_interval != 0 ||
            hpwl_history_.size() < static_cast<std::size_t>(
                std::max(4, config_.adaptive_window / 2))) {
            return;
        }
        const std::size_t half = std::max<std::size_t>(2, hpwl_history_.size() / 2);
        const std::size_t split = hpwl_history_.size() - half;
        const Real old_hpwl = hpwl_history_[split - 1];
        const Real old_overflow = overflow_history_[split - 1];
        const Real hpwl_trend = (hpwl - old_hpwl) /
            std::max<Real>(std::abs(old_hpwl), 1.0);
        const Real overflow_trend = overflow - old_overflow;
        const Real band = std::max<Real>(
            config_.stop_overflow - config_.lower_overflow, 0.005);
        const Real pressure = (overflow - config_.stop_overflow) / band;
        Real signal = 0.0;
        if (overflow > config_.stop_overflow + config_.adaptive_deadband) {
            signal += 0.75 * pressure + 0.25 * std::max<Real>(overflow_trend / band, 0.0);
        } else if (overflow < config_.lower_overflow - config_.adaptive_deadband) {
            signal -= 0.65 * (-pressure) +
                0.20 * std::max<Real>(-overflow_trend / band, 0.0);
        } else if (hpwl_trend > 0.0015) {
            signal -= 0.55 * std::min<Real>(hpwl_trend / 0.01, 1.0);
        } else if (hpwl_trend < -0.0015) {
            signal += 0.20 * std::min<Real>(-hpwl_trend / 0.01, 1.0);
        }
        filtered_signal_ = 0.70 * filtered_signal_ + 0.30 * signal;
        const Real log_step = std::clamp(
            0.07 * config_.adaptive_gain * filtered_signal_,
            -config_.adaptive_max_log_step, config_.adaptive_max_log_step);
        scale_ = std::clamp(scale_ * std::exp(log_step),
                            config_.adaptive_min_scale,
                            config_.adaptive_max_scale);
    }

    Real scale() const noexcept { return scale_; }

private:
    const PlaceConfig& config_;
    Real scale_ = 1.0;
    Real filtered_signal_ = 0.0;
    std::deque<Real> hpwl_history_;
    std::deque<Real> overflow_history_;
};

}  // namespace

PlaceResult global_place(Database& db, const PlaceConfig& config) {
    if (config.iterations <= 0 || config.bins_x <= 0 || config.bins_y <= 0 ||
        config.hpwl_epsilon < 0.0 || config.density_epsilon < 0.0 ||
        config.active_power <= 0.0 || config.density_active_power <= 0.0 ||
        config.step_fraction <= 0.0 ||
        config.stop_overflow < 0.0 || config.stop_overflow >= 1.0) {
        throw std::invalid_argument("invalid global placement configuration");
    }
    if (config.net_batch.weight < 0.0 || config.net_batch.degree_limit < 2 ||
        config.regional_price.update_interval <= 0 ||
        config.regional_price.rho < 0.0 ||
        config.regional_price.min_price <= 0.0 ||
        config.regional_price.max_price < config.regional_price.min_price) {
        throw std::invalid_argument("invalid batch or regional-price configuration");
    }
    if (config.density_net_share.weight < 0.0 ||
        config.density_net_share.degree_limit < 2 ||
        config.density_net_share.zero_gradient_threshold < 0.0 ||
        config.density_net_share.hpwl_orthogonal_projection < 0.0 ||
        config.density_net_share.until_iteration < -1) {
        throw std::invalid_argument("invalid density net sharing configuration");
    }
    if (config.batch_acceptance.max_infeasible_hpwl_increase < 0.0) {
        throw std::invalid_argument("invalid infeasible batch HPWL budget");
    }
    if (config.late_stage.switch_iteration < -1 ||
        config.late_stage.switch_iteration >= config.iterations ||
        config.late_stage.step_multiplier <= 0.0 ||
        config.late_stage.lambda_multiplier <= 0.0) {
        throw std::invalid_argument("invalid late-stage switch configuration");
    }
    configure_threads(config.threads);
    std::filesystem::create_directories(config.output_dir);
    if (config.snapshot_every > 0) {
        std::filesystem::create_directories(config.output_dir / "snapshots");
    }

    const auto started = std::chrono::steady_clock::now();
    PlaceResult result;
    ExactHpwl hpwl_oracle(db);
    ExactOverlapDensity density_oracle(
        db, config.bins_x, config.bins_y, config.target_density);
    result.bisection = recursive_hypergraph_bisection(
        db, density_oracle, config.bisection);
    if (config.bisection.enabled) {
        result.objective_evaluations += 4;
        std::ofstream bisection_metrics(config.output_dir / "bisection_metrics.csv");
        if (!bisection_metrics) {
            throw std::runtime_error("cannot create bisection_metrics.csv");
        }
        bisection_metrics << "splits,leaves,max_depth,fm_moves,cut_net_weight,"
                             "maximum_capacity_error,initial_hpwl,initial_overflow,"
                             "final_hpwl,final_overflow,wall_seconds\n"
                          << std::setprecision(12)
                          << result.bisection.splits << ','
                          << result.bisection.leaves << ','
                          << result.bisection.max_depth << ','
                          << result.bisection.fm_moves << ','
                          << result.bisection.cut_net_weight << ','
                          << result.bisection.maximum_capacity_error << ','
                          << result.bisection.initial_hpwl << ','
                          << result.bisection.initial_overflow << ','
                          << result.bisection.final_hpwl << ','
                          << result.bisection.final_overflow << ','
                          << result.bisection.wall_seconds << '\n';
        std::cout << "[Bisection] splits=" << result.bisection.splits
                  << " leaves=" << result.bisection.leaves
                  << " fm_moves=" << result.bisection.fm_moves
                  << " hpwl=" << result.bisection.initial_hpwl
                  << " -> " << result.bisection.final_hpwl
                  << " overflow=" << result.bisection.initial_overflow
                  << " -> " << result.bisection.final_overflow << '\n';
    }
    result.coarse_flow = coarse_capacity_flow(
        db, density_oracle, config.coarse_flow);
    if (config.coarse_flow.enabled && config.coarse_flow.passes > 0) {
        result.objective_evaluations += 4;
        std::ofstream flow_metrics(config.output_dir / "coarse_flow_metrics.csv");
        if (!flow_metrics) {
            throw std::runtime_error("cannot create coarse_flow_metrics.csv");
        }
        flow_metrics << "passes,flow_edges,moves,exact_anchor_trials,"
                        "exact_anchor_rejections,planned_area,moved_area,"
                        "initial_coarse_overflow,final_coarse_overflow,"
                        "initial_hpwl,initial_overflow,final_hpwl,final_overflow,"
                        "wall_seconds\n"
                     << std::setprecision(12)
                     << result.coarse_flow.passes << ','
                     << result.coarse_flow.flow_edges << ','
                     << result.coarse_flow.moves << ','
                     << result.coarse_flow.exact_anchor_trials << ','
                     << result.coarse_flow.exact_anchor_rejections << ','
                     << result.coarse_flow.planned_area << ','
                     << result.coarse_flow.moved_area << ','
                     << result.coarse_flow.initial_coarse_overflow << ','
                     << result.coarse_flow.final_coarse_overflow << ','
                     << result.coarse_flow.initial_hpwl << ','
                     << result.coarse_flow.initial_overflow << ','
                     << result.coarse_flow.final_hpwl << ','
                     << result.coarse_flow.final_overflow << ','
                     << result.coarse_flow.wall_seconds << '\n';
        std::cout << "[CoarseFlow] passes=" << result.coarse_flow.passes
                  << " edges=" << result.coarse_flow.flow_edges
                  << " moves=" << result.coarse_flow.moves
                  << " exact_anchors=" << result.coarse_flow.exact_anchor_trials
                  << " rejected=" << result.coarse_flow.exact_anchor_rejections
                  << " coarse_overflow="
                  << result.coarse_flow.initial_coarse_overflow
                  << " -> " << result.coarse_flow.final_coarse_overflow
                  << " exact_overflow=" << result.coarse_flow.initial_overflow
                  << " -> " << result.coarse_flow.final_overflow
                  << " hpwl=" << result.coarse_flow.initial_hpwl
                  << " -> " << result.coarse_flow.final_hpwl << '\n';
    }
    result.transport = transport_excess_to_capacity(
        db, density_oracle, config.transport);
    if ((config.transport.rounds > 0 && config.transport.max_moves > 0) ||
        config.transport.identity_exchange_passes > 0) {
        result.objective_evaluations +=
            4 + 2 * config.transport.identity_exchange_passes;
        std::ofstream transport_metrics(config.output_dir / "transport_metrics.csv");
        if (!transport_metrics) {
            throw std::runtime_error("cannot create transport_metrics.csv");
        }
        transport_metrics << "rounds,moves,groups,identity_attempts,"
                             "identity_accepted,identity_permuted_nodes,"
                             "identity_initial_hpwl,identity_final_hpwl,"
                             "identity_initial_overflow,identity_final_overflow,"
                             "identity_wall_seconds,initial_hpwl,initial_overflow,"
                             "final_hpwl,final_overflow,wall_seconds\n"
                          << std::setprecision(12)
                          << result.transport.rounds << ',' << result.transport.moves << ','
                          << result.transport.groups << ','
                          << result.transport.identity_exchange_attempts << ','
                          << result.transport.identity_exchange_accepted << ','
                          << result.transport.identity_exchange_permuted_nodes << ','
                          << result.transport.identity_initial_hpwl << ','
                          << result.transport.identity_final_hpwl << ','
                          << result.transport.identity_initial_overflow << ','
                          << result.transport.identity_final_overflow << ','
                          << result.transport.identity_wall_seconds << ','
                          << result.transport.initial_hpwl << ','
                          << result.transport.initial_overflow << ','
                          << result.transport.final_hpwl << ','
                          << result.transport.final_overflow << ','
                          << result.transport.wall_seconds << '\n';
        std::cout << "[Transport] rounds=" << result.transport.rounds
                  << " moves=" << result.transport.moves
                  << " groups=" << result.transport.groups
                  << " atomic=" << result.transport.atomic_groups
                  << '/' << result.transport.atomic_attempts
                  << " fallback=" << result.transport.atomic_fallbacks
                  << " capacity_reject="
                  << result.transport.atomic_capacity_rejections
                  << " bids=" << result.transport.atomic_feasible_bids
                  << '/' << result.transport.atomic_bid_trials
                  << " rescued=" << result.transport.atomic_rescued_groups
                  << " batch=" << result.transport.atomic_batch_commits
                  << '/' << result.transport.atomic_batch_groups
                  << " identity=" << result.transport.identity_exchange_accepted
                  << '/' << result.transport.identity_exchange_attempts
                  << " identity_hpwl=" << result.transport.identity_initial_hpwl
                  << " -> " << result.transport.identity_final_hpwl
                  << " identity_overflow=" << result.transport.identity_initial_overflow
                  << " -> " << result.transport.identity_final_overflow
                  << " destination="
                  << transport_destination_mode_name(config.transport.destination_mode)
                  << " overflow=" << result.transport.initial_overflow
                  << " -> " << result.transport.final_overflow
                  << " hpwl=" << result.transport.initial_hpwl
                  << " -> " << result.transport.final_hpwl << '\n';
    }
    result.density_coordinate = coordinate_descent_overlap(
        db, density_oracle, config.density_coordinate);
    result.objective_evaluations +=
        result.density_coordinate.sweeps > 0
            ? 2 * result.density_coordinate.sweeps + 2 : 0;
    if (config.density_coordinate.sweeps > 0) {
        std::ofstream coordinate_metrics(
            config.output_dir / "density_coordinate_metrics.csv");
        if (!coordinate_metrics) {
            throw std::runtime_error(
                "cannot create density_coordinate_metrics.csv");
        }
        coordinate_metrics << "sweep,moves,exact_hpwl,overflow\n"
                           << std::setprecision(12);
        for (const DensityCoordinateSweepStats& sweep :
             result.density_coordinate.trajectory) {
            coordinate_metrics << sweep.sweep << ',' << sweep.moves << ','
                               << sweep.hpwl << ',' << sweep.overflow << '\n';
        }
        std::cout << "[DensityCoordinate] sweeps="
                  << result.density_coordinate.sweeps
                  << " moves=" << result.density_coordinate.moves
                  << " hpwl=" << result.density_coordinate.initial_hpwl
                  << " -> " << result.density_coordinate.final_hpwl
                  << " overflow=" << result.density_coordinate.initial_overflow
                  << " -> " << result.density_coordinate.final_overflow
                  << " clusters=" << result.density_coordinate.cluster_moves
                  << '/' << result.density_coordinate.cluster_proposals
                  << " batches="
                  << result.density_coordinate.cooperative_batches
                  << '/' << result.density_coordinate.cooperative_batch_trials
                  << " candidates=" << result.density_coordinate.candidates
                  << '\n';
    }
    result.recovery = recover_hpwl_under_overflow(
        db, density_oracle, config.recovery);
    result.objective_evaluations += result.recovery.objective_evaluations;
    if (config.recovery.sweeps > 0) {
        std::ofstream recovery_metrics(config.output_dir / "recovery_metrics.csv");
        if (!recovery_metrics) {
            throw std::runtime_error("cannot create recovery_metrics.csv");
        }
        recovery_metrics << "sweep,moves,exact_hpwl,overflow\n"
                         << std::setprecision(12);
        for (const RecoverySweepStats& sweep : result.recovery.trajectory) {
            recovery_metrics << sweep.sweep << ',' << sweep.moves << ','
                             << sweep.hpwl << ',' << sweep.overflow << '\n';
        }
        std::cout << "[Recovery] sweeps=" << result.recovery.sweeps
                  << " moves=" << result.recovery.moves
                  << " hpwl=" << result.recovery.initial_hpwl
                  << " -> " << result.recovery.final_hpwl
                  << " overflow=" << result.recovery.initial_overflow
                  << " -> " << result.recovery.final_overflow << '\n';
    }
    result.compact_recovery = compact_support_contraction(
        db, density_oracle, config.compact_recovery);
    result.objective_evaluations += result.compact_recovery.objective_evaluations;
    if (config.compact_recovery.sweeps > 0) {
        std::ofstream compact_metrics(config.output_dir /
                                      "compact_recovery_metrics.csv");
        if (!compact_metrics) {
            throw std::runtime_error("cannot create compact_recovery_metrics.csv");
        }
        compact_metrics << "sweep,moves,exact_hpwl,overflow,compactness\n"
                        << std::setprecision(12);
        for (const CompactRecoverySweepStats& sweep :
             result.compact_recovery.trajectory) {
            compact_metrics << sweep.sweep << ',' << sweep.moves << ','
                            << sweep.hpwl << ',' << sweep.overflow << ','
                            << sweep.compactness << '\n';
        }
        std::cout << "[CompactRecovery] sweeps="
                  << result.compact_recovery.sweeps
                  << " accepted=" << result.compact_recovery.accepted_candidates
                  << " hpwl=" << result.compact_recovery.initial_hpwl
                  << " -> " << result.compact_recovery.final_hpwl
                  << " overflow=" << result.compact_recovery.initial_overflow
                  << " -> " << result.compact_recovery.final_overflow << '\n';
    }
    result.swap_recovery = recover_hpwl_with_equal_shape_swaps(
        db, density_oracle, config.swap_recovery);
    result.objective_evaluations += result.swap_recovery.objective_evaluations;
    if (config.swap_recovery.sweeps > 0 ||
        config.swap_recovery.permutation_sweeps > 0) {
        std::ofstream swap_metrics(config.output_dir / "swap_recovery_metrics.csv");
        if (!swap_metrics) {
            throw std::runtime_error("cannot create swap_recovery_metrics.csv");
        }
        swap_metrics << "sweep,swaps,exact_hpwl,overflow\n"
                     << std::setprecision(12);
        for (const SwapSweepStats& sweep : result.swap_recovery.trajectory) {
            swap_metrics << sweep.sweep << ',' << sweep.swaps << ','
                         << sweep.hpwl << ',' << sweep.overflow << '\n';
        }
        std::cout << "[SwapRecovery] sweeps=" << result.swap_recovery.sweeps
                  << " swaps=" << result.swap_recovery.swaps
                  << " permutation_sweeps="
                  << result.swap_recovery.permutation_sweeps
                  << " hpwl=" << result.swap_recovery.initial_hpwl
                  << " -> " << result.swap_recovery.final_hpwl
                  << " overflow=" << result.swap_recovery.initial_overflow
                  << " -> " << result.swap_recovery.final_overflow << '\n';
    }
    result.post_recovery = recover_hpwl_under_overflow(
        db, density_oracle, config.post_recovery);
    result.objective_evaluations += result.post_recovery.objective_evaluations;
    if (config.post_recovery.sweeps > 0) {
        std::ofstream post_metrics(config.output_dir / "post_recovery_metrics.csv");
        if (!post_metrics) {
            throw std::runtime_error("cannot create post_recovery_metrics.csv");
        }
        post_metrics << "sweep,moves,exact_hpwl,overflow\n"
                     << std::setprecision(12);
        for (const RecoverySweepStats& sweep : result.post_recovery.trajectory) {
            post_metrics << sweep.sweep << ',' << sweep.moves << ','
                         << sweep.hpwl << ',' << sweep.overflow << '\n';
        }
        std::cout << "[PostRecovery] sweeps=" << result.post_recovery.sweeps
                  << " moves=" << result.post_recovery.moves
                  << " hpwl=" << result.post_recovery.initial_hpwl
                  << " -> " << result.post_recovery.final_hpwl
                  << " overflow=" << result.post_recovery.initial_overflow
                  << " -> " << result.post_recovery.final_overflow << '\n';
    }
    result.post_swap_recovery = recover_hpwl_with_equal_shape_swaps(
        db, density_oracle, config.post_swap_recovery);
    result.objective_evaluations += result.post_swap_recovery.objective_evaluations;
    if (config.post_swap_recovery.sweeps > 0) {
        std::ofstream post_swap_metrics(
            config.output_dir / "post_swap_recovery_metrics.csv");
        if (!post_swap_metrics) {
            throw std::runtime_error("cannot create post_swap_recovery_metrics.csv");
        }
        post_swap_metrics << "sweep,swaps,exact_hpwl,overflow\n"
                          << std::setprecision(12);
        for (const SwapSweepStats& sweep : result.post_swap_recovery.trajectory) {
            post_swap_metrics << sweep.sweep << ',' << sweep.swaps << ','
                              << sweep.hpwl << ',' << sweep.overflow << '\n';
        }
        std::cout << "[PostSwapRecovery] sweeps="
                  << result.post_swap_recovery.sweeps
                  << " swaps=" << result.post_swap_recovery.swaps
                  << " hpwl=" << result.post_swap_recovery.initial_hpwl
                  << " -> " << result.post_swap_recovery.final_hpwl
                  << " overflow=" << result.post_swap_recovery.initial_overflow
                  << " -> " << result.post_swap_recovery.final_overflow << '\n';
    }
    result.assignment_recovery = recover_hpwl_with_equal_shape_swaps(
        db, density_oracle, config.assignment_recovery);
    result.objective_evaluations += result.assignment_recovery.objective_evaluations;
    if (config.assignment_recovery.assignment_sweeps > 0) {
        std::ofstream assignment_metrics(
            config.output_dir / "anchor_assignment_metrics.csv");
        if (!assignment_metrics) {
            throw std::runtime_error("cannot create anchor_assignment_metrics.csv");
        }
        assignment_metrics << "sweep,moved_nodes,exact_hpwl,overflow\n"
                           << std::setprecision(12);
        for (const SwapSweepStats& sweep : result.assignment_recovery.trajectory) {
            assignment_metrics << sweep.sweep << ',' << sweep.swaps << ','
                               << sweep.hpwl << ',' << sweep.overflow << '\n';
        }
        std::cout << "[AnchorAssignment] sweeps="
                  << result.assignment_recovery.assignment_sweeps
                  << " batches=" << result.assignment_recovery.assignment_batches
                  << " moved_nodes=" << result.assignment_recovery.assigned_nodes
                  << " hpwl=" << result.assignment_recovery.initial_hpwl
                  << " -> " << result.assignment_recovery.final_hpwl
                  << " overflow=" << result.assignment_recovery.initial_overflow
                  << " -> " << result.assignment_recovery.final_overflow << '\n';
    }
    // Overlap is an exact rectangle-bin objective.  A zero configuration
    // value must remain zero; using a bin pitch here silently changes the
    // objective and makes the logged overflow differ from the audited one.
    const Real base_density_epsilon = config.density_epsilon;
    LambdaConfig lambda_config = config.lambda;
    lambda_config.stop_overflow = config.stop_overflow;
    LambdaController lambda(lambda_config);
    std::unique_ptr<Optimizer> optimizer = make_optimizer(
        config.optimizer, config.beta1, config.beta2,
        config.momentum, config.numerical_epsilon);
    optimizer->reset(2 * db.movable_ids.size());
    AdaptiveEpsilon adaptive(config);

    std::ofstream metrics(config.output_dir / "global_metrics.csv");
    if (!metrics) throw std::runtime_error("cannot create global_metrics.csv");
    metrics << "iteration,exact_hpwl,overflow,max_density,density_energy,"
               "lambda_base,lambda_control,lambda_effective,hpwl_epsilon,"
               "density_epsilon,epsilon_scale,learning_rate,optimizer,"
               "batch_accepted,batch_scale,batch_trials,batch_rejected\n";
    metrics << std::setprecision(12);

    const std::size_t n = db.movable_ids.size();
    const Real base_learning_rate = config.step_fraction *
        std::min(db.xh - db.xl, db.yh - db.yl);
    std::vector<Real> wire_x, wire_y, density_x, density_y;
    std::vector<Real> batch_x, batch_y;
    std::vector<Real> regional_prices = config.regional_price.initial_prices;
    const std::size_t regional_size =
        static_cast<std::size_t>(config.bins_x) * config.bins_y;
    if (regional_prices.empty()) {
        regional_prices.assign(regional_size, 1.0);
    } else if (regional_prices.size() != regional_size) {
        throw std::invalid_argument(
            "initial regional price field has wrong dimensions");
    }
    std::vector<Real> gradient(2 * n), delta;
    std::vector<Real> positions;
    std::vector<Real> best_feasible;
    std::vector<Real> best_feasible_prices;
    std::vector<Real> best_overflow = capture(db);
    std::vector<Real> best_overflow_prices = regional_prices;
    Real best_feasible_hpwl = std::numeric_limits<Real>::infinity();
    Real best_overflow_value = std::numeric_limits<Real>::infinity();
    int best_feasible_iteration = -1;
    // Recovery and transport are exact stages that precede GP. Preserve an
    // exact-feasible checkpoint so a later unconstrained GP step cannot
    // overwrite it with a lower-HPWL but infeasible state.
    const std::vector<Real> pre_gp_positions = capture(db);
    const DensityMetrics pre_gp_density = density_oracle.evaluate(
        0.0, config.density_active_power, nullptr, nullptr);
    const Real pre_gp_hpwl = hpwl_oracle.evaluate(
        0.0, config.active_power, -1, nullptr, nullptr);
    if (pre_gp_density.overflow < best_overflow_value) {
        best_overflow_value = pre_gp_density.overflow;
        best_overflow = pre_gp_positions;
        best_overflow_prices = regional_prices;
    }
    if (pre_gp_density.overflow <= config.stop_overflow) {
        best_feasible_hpwl = pre_gp_hpwl;
        best_feasible = pre_gp_positions;
        best_feasible_prices = regional_prices;
    }
    bool lambda_initialized = false;
    IterationMetrics last_metrics;

    for (int iteration = 0; iteration < config.iterations; ++iteration) {
        const bool late_stage = config.late_stage.switch_iteration >= 0 &&
            iteration >= config.late_stage.switch_iteration;
        const Real step_multiplier = late_stage
            ? config.late_stage.step_multiplier : 1.0;
        const Real lambda_multiplier = late_stage
            ? config.late_stage.lambda_multiplier : 1.0;
        const Real learning_rate = base_learning_rate * step_multiplier;
        if (late_stage && iteration == config.late_stage.switch_iteration &&
            config.late_stage.reset_optimizer) {
            optimizer->reset(2 * db.movable_ids.size());
        }
        const Real epsilon_scale = adaptive.scale();
        const Real hpwl_epsilon = config.hpwl_epsilon * epsilon_scale;
        const Real density_epsilon = base_density_epsilon * epsilon_scale;
        const DensityMetrics density = config.regional_price.enabled
            ? density_oracle.evaluate_with_prices(
                regional_prices, density_epsilon, config.density_active_power,
                &density_x, &density_y)
            : density_oracle.evaluate(
                density_epsilon, config.density_active_power, &density_x, &density_y);
        DensityNetShareConfig sharing = config.density_net_share;
        if (sharing.until_iteration >= 0 && iteration >= sharing.until_iteration) {
            sharing.weight = 0.0;
            sharing.hpwl_orthogonal_projection = 0.0;
        }
        share_density_direction_over_nets(db, sharing, density_x, density_y);
        (void)hpwl_oracle.evaluate(
            hpwl_epsilon, config.active_power, config.degree_limit,
            &wire_x, &wire_y);
        project_density_direction_from_hpwl(
            db, sharing.hpwl_orthogonal_projection,
            wire_x, wire_y, density_x, density_y);
        const DensityMetrics exact_density = density_epsilon == 0.0
            ? density
            : density_oracle.evaluate(0.0, config.density_active_power,
                                      nullptr, nullptr);
        const Real exact_hpwl = hpwl_oracle.evaluate(
            0.0, 1.0, -1, nullptr, nullptr);
        result.objective_evaluations += 2;
        if (!lambda_initialized) {
            lambda.initialize(db, wire_x, wire_y, density_x, density_y,
                              exact_hpwl, exact_density.overflow);
            lambda_initialized = true;
        } else {
            lambda.update(iteration, exact_hpwl, exact_density.overflow);
        }

        last_metrics.iteration = iteration;
        last_metrics.hpwl = exact_hpwl;
        last_metrics.overflow = exact_density.overflow;
        last_metrics.max_density = exact_density.max_density;
        last_metrics.density_energy = exact_density.energy;
        last_metrics.lambda_base = lambda.base();
        last_metrics.lambda_control = lambda.control();
        last_metrics.lambda_effective = lambda.effective() * lambda_multiplier;
        last_metrics.hpwl_epsilon = hpwl_epsilon;
        last_metrics.density_epsilon = density_epsilon;
        last_metrics.epsilon_scale = epsilon_scale;
        last_metrics.learning_rate = learning_rate;
        positions = capture(db);
        if (exact_density.overflow < best_overflow_value) {
            best_overflow_value = exact_density.overflow;
            best_overflow = positions;
            best_overflow_prices = regional_prices;
        }
        if (exact_density.overflow <= config.stop_overflow &&
            exact_hpwl < best_feasible_hpwl) {
            best_feasible_hpwl = exact_hpwl;
            best_feasible = positions;
            best_feasible_prices = regional_prices;
            best_feasible_iteration = iteration;
        }
        save_snapshot(db, config, iteration);
        if (iteration % config.log_every == 0 || iteration + 1 == config.iterations) {
            std::cout << "[GP] iter=" << iteration << " hpwl=" << exact_hpwl
                      << " overflow=" << exact_density.overflow
                      << " lambda=" << lambda.effective() * lambda_multiplier
                      << " epsilon_scale=" << epsilon_scale << '\n';
        }

        const Real effective_lambda = lambda.effective() * lambda_multiplier;
        if (config.regional_price.enabled &&
            (iteration + 1) % config.regional_price.update_interval == 0) {
            update_regional_prices(density_oracle, config.regional_price,
                                   regional_prices);
        }
        if (config.net_batch.weight > 0.0) {
            compute_net_batch_direction(
                db, density_x, density_y, config.net_batch.degree_limit,
                batch_x, batch_y);
        } else {
            batch_x.assign(db.nodes.size(), 0.0);
            batch_y.assign(db.nodes.size(), 0.0);
        }
        #pragma omp parallel for schedule(static)
        for (int movable = 0; movable < static_cast<int>(n); ++movable) {
            const int id = db.movable_ids[movable];
            const Node& node = db.nodes[id];
            const Real preconditioner = std::max<Real>(
                1.0, db.node_pin_count[id] + effective_lambda * node.area());
            gradient[movable] =
                (wire_x[id] + effective_lambda *
                    (density_x[id] + config.net_batch.weight * batch_x[id])) /
                preconditioner;
            gradient[n + movable] =
                (wire_y[id] + effective_lambda *
                    (density_y[id] + config.net_batch.weight * batch_y[id])) /
                preconditioner;
        }
        BatchAcceptanceResult batch;
        if (iteration + 1 < config.iterations) {
            optimizer->compute_delta(
                gradient, learning_rate,
                config.max_step_multiplier * learning_rate, delta);
            if (config.batch_acceptance.enabled) {
                BatchAcceptanceConfig acceptance = config.batch_acceptance;
                acceptance.overflow_cap = config.stop_overflow;
                batch = accept_exact_batch(
                    positions, delta, exact_hpwl, exact_density.overflow, acceptance,
                    [&](const std::vector<Real>& candidate) {
                        apply_positions(db, candidate);
                    },
                    [&]() {
                        const Real exact_hpwl = hpwl_oracle.evaluate(
                            0.0, config.active_power, -1, nullptr, nullptr);
                        const DensityMetrics exact_density = density_oracle.evaluate(
                            0.0, config.density_active_power, nullptr, nullptr);
                        return std::make_pair(exact_hpwl, exact_density.overflow);
                    });
                result.objective_evaluations += 2 * batch.trials;
                result.batch_trials += batch.trials;
                if (batch.accepted) {
                    ++result.accepted_batches;
                } else {
                    ++result.rejected_batches;
                }
            } else {
                #pragma omp parallel for schedule(static)
                for (int coordinate = 0;
                     coordinate < static_cast<int>(positions.size()); ++coordinate) {
                    positions[coordinate] -= delta[coordinate];
                }
                apply_positions(db, positions);
                batch.accepted = true;
                batch.scale = 1.0;
            }
        }
        metrics << iteration << ',' << exact_hpwl << ',' << exact_density.overflow << ','
                << exact_density.max_density << ',' << exact_density.energy << ','
                << lambda.base() << ',' << lambda.control() << ','
                << lambda.effective() * lambda_multiplier << ',' << hpwl_epsilon << ','
                << density_epsilon << ',' << epsilon_scale << ','
                << learning_rate << ',' << optimizer->name() << ','
                << (batch.accepted ? 1 : 0) << ',' << batch.scale << ','
                << batch.trials << ','
                << ((!batch.accepted && batch.trials > 0) ? 1 : 0) << '\n';
        adaptive.observe(iteration, exact_hpwl, exact_density.overflow);
    }

    // A staged experiment may continue only from a checkpoint produced by
    // this run.  Keep the last trajectory state separate from the feasible
    // selector below, whose purpose is final reporting rather than handoff.
    write_bookshelf_placement(db, config.output_dir / "last.pl");

    if (!best_feasible.empty()) {
        apply_positions(db, best_feasible);
        result.final_regional_prices = best_feasible_prices;
        result.feasible = true;
        result.selected_iteration = best_feasible_iteration;
    } else {
        apply_positions(db, best_overflow);
        result.final_regional_prices = best_overflow_prices;
        result.feasible = false;
        result.selected_iteration = -1;
    }
    const DensityMetrics final_density = density_oracle.evaluate(
        0.0, config.density_active_power, nullptr, nullptr);
    const Real final_hpwl = hpwl_oracle.evaluate(
        0.0, config.active_power, -1, nullptr, nullptr);
    result.objective_evaluations += 2;
    result.final_metrics = last_metrics;
    result.final_metrics.iteration = result.selected_iteration;
    result.final_metrics.hpwl = final_hpwl;
    result.final_metrics.overflow = final_density.overflow;
    result.final_metrics.max_density = final_density.max_density;
    result.final_metrics.density_energy = final_density.energy;
    result.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return result;
}

}  // namespace ea
