#include "epsilon_active/recovery.hpp"

#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ea {
namespace {

std::vector<std::vector<int>> build_node_nets(const Database& db) {
    std::vector<std::vector<int>> result(db.nodes.size());
    for (int net_id = 0; net_id < static_cast<int>(db.nets.size()); ++net_id) {
        const Net& net = db.nets[net_id];
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            result[db.pins[p].node].push_back(net_id);
        }
    }
    for (auto& nets : result) {
        std::sort(nets.begin(), nets.end());
        nets.erase(std::unique(nets.begin(), nets.end()), nets.end());
    }
    return result;
}

Real affected_hpwl_delta(const Database& db, int node_id, Real new_x, Real new_y,
                         const std::vector<int>& affected_nets) {
    Real delta = 0.0;
    for (int net_id : affected_nets) {
        const Net& net = db.nets[net_id];
        if (net.pin_count < 2) continue;
        Real old_min_x = std::numeric_limits<Real>::infinity();
        Real old_max_x = -std::numeric_limits<Real>::infinity();
        Real old_min_y = std::numeric_limits<Real>::infinity();
        Real old_max_y = -std::numeric_limits<Real>::infinity();
        Real new_min_x = old_min_x;
        Real new_max_x = old_max_x;
        Real new_min_y = old_min_y;
        Real new_max_y = old_max_y;
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            const Pin& pin = db.pins[p];
            const Node& node = db.nodes[pin.node];
            const Real old_pin_x = node.x + pin.offset_x;
            const Real old_pin_y = node.y + pin.offset_y;
            const Real new_pin_x = (pin.node == node_id ? new_x : node.x) + pin.offset_x;
            const Real new_pin_y = (pin.node == node_id ? new_y : node.y) + pin.offset_y;
            old_min_x = std::min(old_min_x, old_pin_x);
            old_max_x = std::max(old_max_x, old_pin_x);
            old_min_y = std::min(old_min_y, old_pin_y);
            old_max_y = std::max(old_max_y, old_pin_y);
            new_min_x = std::min(new_min_x, new_pin_x);
            new_max_x = std::max(new_max_x, new_pin_x);
            new_min_y = std::min(new_min_y, new_pin_y);
            new_max_y = std::max(new_max_y, new_pin_y);
        }
        delta += net.weight * ((new_max_x - new_min_x) + (new_max_y - new_min_y) -
                               (old_max_x - old_min_x) - (old_max_y - old_min_y));
    }
    return delta;
}

Real affected_hpwl_delta_group(
    const Database& db, const std::vector<int>& node_ids, Real offset_x,
    Real offset_y, const std::vector<int>& affected_nets,
    std::vector<int>& node_marks, int mark) {
    for (int id : node_ids) node_marks[id] = mark;
    Real delta = 0.0;
    for (int net_id : affected_nets) {
        const Net& net = db.nets[net_id];
        if (net.pin_count < 2) continue;
        Real old_min_x = std::numeric_limits<Real>::infinity();
        Real old_max_x = -std::numeric_limits<Real>::infinity();
        Real old_min_y = std::numeric_limits<Real>::infinity();
        Real old_max_y = -std::numeric_limits<Real>::infinity();
        Real new_min_x = old_min_x;
        Real new_max_x = old_max_x;
        Real new_min_y = old_min_y;
        Real new_max_y = old_max_y;
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            const Pin& pin = db.pins[p];
            const Node& node = db.nodes[pin.node];
            const Real old_pin_x = node.x + pin.offset_x;
            const Real old_pin_y = node.y + pin.offset_y;
            const bool moved = node_marks[pin.node] == mark;
            const Real new_pin_x = old_pin_x + (moved ? offset_x : 0.0);
            const Real new_pin_y = old_pin_y + (moved ? offset_y : 0.0);
            old_min_x = std::min(old_min_x, old_pin_x);
            old_max_x = std::max(old_max_x, old_pin_x);
            old_min_y = std::min(old_min_y, old_pin_y);
            old_max_y = std::max(old_max_y, old_pin_y);
            new_min_x = std::min(new_min_x, new_pin_x);
            new_max_x = std::max(new_max_x, new_pin_x);
            new_min_y = std::min(new_min_y, new_pin_y);
            new_max_y = std::max(new_max_y, new_pin_y);
        }
        delta += net.weight *
            ((new_max_x - new_min_x) + (new_max_y - new_min_y) -
             (old_max_x - old_min_x) - (old_max_y - old_min_y));
    }
    return delta;
}

Real affected_hpwl_delta_group_positions(
    const Database& db, const std::vector<DensityNodeMove>& moves,
    const std::vector<int>& affected_nets, std::vector<int>& node_marks,
    std::vector<Real>& trial_x, std::vector<Real>& trial_y, int mark) {
    for (const DensityNodeMove& move : moves) {
        node_marks[move.node_id] = mark;
        trial_x[move.node_id] = move.x;
        trial_y[move.node_id] = move.y;
    }
    Real delta = 0.0;
    for (int net_id : affected_nets) {
        const Net& net = db.nets[net_id];
        if (net.pin_count < 2) continue;
        Real old_min_x = std::numeric_limits<Real>::infinity();
        Real old_max_x = -std::numeric_limits<Real>::infinity();
        Real old_min_y = std::numeric_limits<Real>::infinity();
        Real old_max_y = -std::numeric_limits<Real>::infinity();
        Real new_min_x = old_min_x;
        Real new_max_x = old_max_x;
        Real new_min_y = old_min_y;
        Real new_max_y = old_max_y;
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            const Pin& pin = db.pins[p];
            const Node& node = db.nodes[pin.node];
            const Real old_x = node.x + pin.offset_x;
            const Real old_y = node.y + pin.offset_y;
            const bool moved = node_marks[pin.node] == mark;
            const Real new_x = (moved ? trial_x[pin.node] : node.x) + pin.offset_x;
            const Real new_y = (moved ? trial_y[pin.node] : node.y) + pin.offset_y;
            old_min_x = std::min(old_min_x, old_x);
            old_max_x = std::max(old_max_x, old_x);
            old_min_y = std::min(old_min_y, old_y);
            old_max_y = std::max(old_max_y, old_y);
            new_min_x = std::min(new_min_x, new_x);
            new_max_x = std::max(new_max_x, new_x);
            new_min_y = std::min(new_min_y, new_y);
            new_max_y = std::max(new_max_y, new_y);
        }
        delta += net.weight *
            ((new_max_x - new_min_x) + (new_max_y - new_min_y) -
             (old_max_x - old_min_x) - (old_max_y - old_min_y));
    }
    return delta;
}

std::vector<Real> breakpoint_positions(Real center, Real extent, Real lower,
                                       Real upper, Real bin_width, Real direction) {
    std::vector<Real> positions;
    if (direction == 0.0 || bin_width <= 0.0) return positions;
    const Real edge_a = center - 0.5 * extent;
    const Real edge_b = center + 0.5 * extent;
    const Real first_edge = direction > 0.0 ? edge_a : edge_b;
    const Real edge_limit = direction > 0.0 ? upper + 0.5 * extent
                                            : lower - 0.5 * extent;
    const Real shifted = first_edge - lower;
    const Real base = direction > 0.0
        ? std::floor(shifted / bin_width) + 1.0
        : std::ceil(shifted / bin_width) - 1.0;
    for (int k = 0; k < 4; ++k) {
        const Real boundary = lower + (base + direction * static_cast<Real>(k)) * bin_width;
        const Real distance = direction * (boundary - first_edge);
        if (distance <= 1.0e-8 || direction * (boundary - first_edge) >=
            direction * (edge_limit - first_edge)) continue;
        positions.push_back(std::clamp(center + direction * distance,
                                       lower + 0.5 * extent,
                                       upper - 0.5 * extent));
    }
    return positions;
}

bool boundary_strip_active(Real center, Real extent, Real lower,
                           Real upper, Real bin_width, Real epsilon_bins) {
    if (bin_width <= 0.0 || epsilon_bins < 0.0) return false;
    const Real lo = center - 0.5 * extent;
    const Real hi = center + 0.5 * extent;
    const Real eps = epsilon_bins * bin_width;
    const auto distance = [&](Real edge) {
        const Real q = (edge - lower) / bin_width;
        const Real nearest = lower + std::round(q) * bin_width;
        return std::abs(edge - nearest);
    };
    return distance(lo) <= eps + 1.0e-10 ||
           distance(hi) <= eps + 1.0e-10 ||
           lo <= lower + eps + 1.0e-10 || hi >= upper - eps - 1.0e-10;
}

std::uint64_t mix_noise(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

}  // namespace

RecoveryStats recover_hpwl_under_overflow(
    Database& db, ExactOverlapDensity& density, const RecoveryConfig& config) {
    if (config.sweeps < 0 || config.line_search_steps <= 0 ||
        config.degree_limit < 2 || config.step_bins <= 0.0 ||
        config.hpwl_epsilon < 0.0 || config.active_power <= 0.0 ||
        config.net_batch_weight < 0.0 || config.net_batch_degree_limit < 2 ||
        config.net_block_degree_limit < 2 || config.net_block_max_nodes < 2 ||
        config.net_block_max_blocks < 0 ||
        config.objective_density_weight < 0.0 ||
        config.boundary_epsilon_bins < 0.0) {
        throw std::invalid_argument("invalid constrained recovery configuration");
    }
    RecoveryStats stats;
    if (config.sweeps == 0) return stats;

    const auto started = std::chrono::steady_clock::now();
    ExactHpwl hpwl(db);
    DensityMetrics density_metrics = density.evaluate(0.0, 1.0, nullptr, nullptr);
    stats.initial_overflow = density_metrics.overflow;
    stats.initial_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    stats.objective_evaluations += 2;
    const auto node_nets = build_node_nets(db);
    const Real overflow_tolerance = 1.0e-10 * density.bin_area();
    const bool absolute_overflow_cap = config.overflow_cap >= 0.0;
    if (config.overflow_cap >= 1.0) {
        throw std::invalid_argument("recovery overflow cap must be below one");
    }
    const Real overflow_area_cap = absolute_overflow_cap
        ? config.overflow_cap * db.movable_area
        : stats.initial_overflow * db.movable_area;
    const Real movable_area_scale = std::max<Real>(db.movable_area, 1.0);
    const auto combined_delta = [&](Real hpwl_delta, Real overflow_area_delta) {
        return hpwl_delta + config.objective_density_weight *
            (overflow_area_delta / movable_area_scale);
    };
    Real current_overflow_area = stats.initial_overflow * db.movable_area;
    std::vector<Real> grad_x, grad_y;
    std::vector<Real> exact_density_grad_x, exact_density_grad_y;
    std::vector<Real> batch_x(db.nodes.size(), 0.0);
    std::vector<Real> batch_y(db.nodes.size(), 0.0);
    std::vector<Real> compact_x(db.nodes.size(), 0.0);
    std::vector<Real> compact_y(db.nodes.size(), 0.0);
    std::vector<int> order = db.movable_ids;
    std::vector<int> group_marks(db.nodes.size(), 0);
    std::vector<Real> group_trial_x(db.nodes.size(), 0.0);
    std::vector<Real> group_trial_y(db.nodes.size(), 0.0);
    const std::vector<int> no_nodes;
    int group_mark = 0;

    for (int sweep = 0; sweep < config.sweeps; ++sweep) {
        hpwl.evaluate(config.hpwl_epsilon, config.active_power,
                      config.degree_limit, &grad_x, &grad_y);
        ++stats.objective_evaluations;
        if (config.net_block && config.net_block_density_direction) {
            density.evaluate(0.0, 1.0, &exact_density_grad_x,
                             &exact_density_grad_y);
            ++stats.objective_evaluations;
        }
        if (config.net_batch_weight > 0.0) {
            std::vector<int> counts(db.nodes.size(), 0);
            std::fill(batch_x.begin(), batch_x.end(), 0.0);
            std::fill(batch_y.begin(), batch_y.end(), 0.0);
            for (const Net& net : db.nets) {
                if (net.pin_count < 2 ||
                    static_cast<int>(net.pin_count) > config.net_batch_degree_limit) {
                    continue;
                }
                Real mean_x = 0.0;
                Real mean_y = 0.0;
                int count = 0;
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const int id = db.pins[p].node;
                    if (db.nodes[id].fixed) continue;
                    mean_x += grad_x[id];
                    mean_y += grad_y[id];
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
            for (int id : db.movable_ids) {
                if (counts[id] > 0) {
                    batch_x[id] /= static_cast<Real>(counts[id]);
                    batch_y[id] /= static_cast<Real>(counts[id]);
                }
            }
        }
        Real support_x = 0.0;
        Real support_y = 0.0;
        Real support_weight = 0.0;
        if (config.compact_directions) {
            std::vector<int> counts(db.nodes.size(), 0);
            std::fill(compact_x.begin(), compact_x.end(), 0.0);
            std::fill(compact_y.begin(), compact_y.end(), 0.0);
            for (int id : db.movable_ids) {
                const Node& node = db.nodes[id];
                const Real area = std::max<Real>(node.width * node.height, 1.0);
                support_x += area * node.x;
                support_y += area * node.y;
                support_weight += area;
            }
            support_x /= std::max<Real>(support_weight, 1.0);
            support_y /= std::max<Real>(support_weight, 1.0);
            for (const Net& net : db.nets) {
                if (net.pin_count < 2 ||
                    static_cast<int>(net.pin_count) > config.degree_limit) {
                    continue;
                }
                Real mean_x = 0.0;
                Real mean_y = 0.0;
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const Pin& pin = db.pins[p];
                    mean_x += db.nodes[pin.node].x + pin.offset_x;
                    mean_y += db.nodes[pin.node].y + pin.offset_y;
                }
                mean_x /= static_cast<Real>(net.pin_count);
                mean_y /= static_cast<Real>(net.pin_count);
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const int id = db.pins[p].node;
                    if (db.nodes[id].fixed) continue;
                    compact_x[id] += mean_x;
                    compact_y[id] += mean_y;
                    ++counts[id];
                }
            }
            for (int id : db.movable_ids) {
                if (counts[id] > 0) {
                    compact_x[id] /= static_cast<Real>(counts[id]);
                    compact_y[id] /= static_cast<Real>(counts[id]);
                } else {
                    compact_x[id] = support_x;
                    compact_y[id] = support_y;
                }
            }
        }
        int sweep_block_moves = 0;
        if (config.net_block && config.net_block_max_blocks > 0) {
            struct BlockCandidate {
                std::vector<int> nodes;
                std::vector<int> nets;
                Real direction_x = 0.0;
                Real direction_y = 0.0;
                Real density_direction_x = 0.0;
                Real density_direction_y = 0.0;
                Real target_x = 0.0;
                Real target_y = 0.0;
                Real score = 0.0;
            };
            std::vector<BlockCandidate> blocks;
            blocks.reserve(db.nets.size());
            for (const Net& net : db.nets) {
                if (net.pin_count < 2 ||
                    static_cast<int>(net.pin_count) >
                        config.net_block_degree_limit) {
                    continue;
                }
                BlockCandidate block;
                Real min_pin_x = std::numeric_limits<Real>::infinity();
                Real max_pin_x = -std::numeric_limits<Real>::infinity();
                Real min_pin_y = std::numeric_limits<Real>::infinity();
                Real max_pin_y = -std::numeric_limits<Real>::infinity();
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const Pin& pin = db.pins[p];
                    const int id = pin.node;
                    const Real pin_x = db.nodes[id].x + pin.offset_x;
                    const Real pin_y = db.nodes[id].y + pin.offset_y;
                    min_pin_x = std::min(min_pin_x, pin_x);
                    max_pin_x = std::max(max_pin_x, pin_x);
                    min_pin_y = std::min(min_pin_y, pin_y);
                    max_pin_y = std::max(max_pin_y, pin_y);
                    if (!db.nodes[id].fixed) block.nodes.push_back(id);
                }
                block.target_x = 0.5 * (min_pin_x + max_pin_x);
                block.target_y = 0.5 * (min_pin_y + max_pin_y);
                std::sort(block.nodes.begin(), block.nodes.end());
                block.nodes.erase(
                    std::unique(block.nodes.begin(), block.nodes.end()),
                    block.nodes.end());
                if (block.nodes.size() < 2 ||
                    block.nodes.size() >
                        static_cast<std::size_t>(config.net_block_max_nodes)) {
                    continue;
                }
                for (int id : block.nodes) {
                    block.direction_x += grad_x[id];
                    block.direction_y += grad_y[id];
                    if (config.net_block_density_direction) {
                        block.density_direction_x += exact_density_grad_x[id];
                        block.density_direction_y += exact_density_grad_y[id];
                    }
                    block.nets.insert(block.nets.end(), node_nets[id].begin(),
                                      node_nets[id].end());
                }
                std::sort(block.nets.begin(), block.nets.end());
                block.nets.erase(std::unique(block.nets.begin(), block.nets.end()),
                                 block.nets.end());
                block.score =
                    (std::hypot(block.direction_x, block.direction_y) +
                     std::hypot(block.density_direction_x,
                                block.density_direction_y)) /
                    std::sqrt(static_cast<Real>(block.nodes.size()));
                if (config.net_block_contraction) {
                    block.score += (max_pin_x - min_pin_x) +
                                   (max_pin_y - min_pin_y);
                }
                if (block.score > 0.0) blocks.push_back(std::move(block));
            }
            std::sort(blocks.begin(), blocks.end(),
                      [](const BlockCandidate& a, const BlockCandidate& b) {
                          return a.score > b.score;
                      });
            std::vector<unsigned char> used(db.nodes.size(), 0);
            int processed_blocks = 0;
            for (const BlockCandidate& block : blocks) {
                if (processed_blocks >= config.net_block_max_blocks) break;
                bool disjoint = true;
                for (int id : block.nodes) {
                    if (used[id]) {
                        disjoint = false;
                        break;
                    }
                }
                if (!disjoint) continue;
                for (int id : block.nodes) used[id] = 1;
                ++processed_blocks;
                ++stats.net_block_attempts;

                Real min_offset_x = -std::numeric_limits<Real>::infinity();
                Real max_offset_x = std::numeric_limits<Real>::infinity();
                Real min_offset_y = -std::numeric_limits<Real>::infinity();
                Real max_offset_y = std::numeric_limits<Real>::infinity();
                for (int id : block.nodes) {
                    const Node& node = db.nodes[id];
                    min_offset_x = std::max(
                        min_offset_x, db.xl + 0.5 * node.width - node.x);
                    max_offset_x = std::min(
                        max_offset_x, db.xh - 0.5 * node.width - node.x);
                    min_offset_y = std::max(
                        min_offset_y, db.yl + 0.5 * node.height - node.y);
                    max_offset_y = std::min(
                        max_offset_y, db.yh - 0.5 * node.height - node.y);
                }

                bool found = false;
                Real best_offset_x = 0.0;
                Real best_offset_y = 0.0;
                Real best_hpwl_delta = 0.0;
                Real best_objective_delta = 0.0;
                DensityMove best_density_move;
                std::vector<DensityNodeMove> best_contraction_moves;
                auto consider_block = [&](Real offset_x, Real offset_y) {
                    offset_x = std::clamp(offset_x, min_offset_x, max_offset_x);
                    offset_y = std::clamp(offset_y, min_offset_y, max_offset_y);
                    if (std::abs(offset_x) <= 1.0e-12 &&
                        std::abs(offset_y) <= 1.0e-12) {
                        return;
                    }
                    ++group_mark;
                    const Real hpwl_delta = affected_hpwl_delta_group(
                        db, block.nodes, offset_x, offset_y, block.nets,
                        group_marks, group_mark);
                    std::vector<DensityNodeMove> moves;
                    moves.reserve(block.nodes.size());
                    for (int id : block.nodes) {
                        moves.push_back(
                            {id, db.nodes[id].x + offset_x,
                             db.nodes[id].y + offset_y});
                    }
                    DensityMove density_move =
                        density.evaluate_group_move(moves);
                    const Real objective_delta = combined_delta(
                        hpwl_delta, density_move.overflow_area_delta);
                    if (objective_delta >= -1.0e-12) return;
                    if (absolute_overflow_cap) {
                        if (current_overflow_area >
                            overflow_area_cap + overflow_tolerance) {
                            if (density_move.overflow_area_delta >
                                overflow_tolerance) {
                                return;
                            }
                        } else if (current_overflow_area +
                                   density_move.overflow_area_delta >
                                   overflow_area_cap + overflow_tolerance) {
                            return;
                        }
                    } else if (density_move.overflow_area_delta >
                               overflow_tolerance) {
                        return;
                    }
                    if (!found || objective_delta < best_objective_delta) {
                        found = true;
                        best_offset_x = offset_x;
                        best_offset_y = offset_y;
                        best_hpwl_delta = hpwl_delta;
                        best_objective_delta = objective_delta;
                        best_density_move = std::move(density_move);
                        best_contraction_moves.clear();
                    }
                };
                auto consider_contraction = [&](Real alpha, int mode) {
                    std::vector<DensityNodeMove> moves;
                    moves.reserve(block.nodes.size());
                    for (int id : block.nodes) {
                        const Node& node = db.nodes[id];
                        Real new_x = node.x + alpha * (block.target_x - node.x);
                        Real new_y = node.y + alpha * (block.target_y - node.y);
                        if (mode == 1) new_y = node.y;
                        if (mode == 2) new_x = node.x;
                        new_x = std::clamp(new_x,
                            db.xl + 0.5 * node.width,
                            db.xh - 0.5 * node.width);
                        new_y = std::clamp(new_y,
                            db.yl + 0.5 * node.height,
                            db.yh - 0.5 * node.height);
                        moves.push_back({id, new_x, new_y});
                    }
                    ++stats.net_block_contraction_attempts;
                    ++group_mark;
                    const Real hpwl_delta = affected_hpwl_delta_group_positions(
                        db, moves, block.nets, group_marks, group_trial_x,
                        group_trial_y, group_mark);
                    DensityMove density_move = density.evaluate_group_move(moves);
                    const Real objective_delta = combined_delta(
                        hpwl_delta, density_move.overflow_area_delta);
                    if (objective_delta >= -1.0e-12) return;
                    if (absolute_overflow_cap) {
                        if (current_overflow_area >
                            overflow_area_cap + overflow_tolerance) {
                            if (density_move.overflow_area_delta >
                                overflow_tolerance) return;
                        } else if (current_overflow_area +
                                   density_move.overflow_area_delta >
                                   overflow_area_cap + overflow_tolerance) {
                            return;
                        }
                    } else if (density_move.overflow_area_delta >
                               overflow_tolerance) {
                        return;
                    }
                    if (!found || objective_delta < best_objective_delta) {
                        found = true;
                        best_hpwl_delta = hpwl_delta;
                        best_objective_delta = objective_delta;
                        best_density_move = std::move(density_move);
                        best_contraction_moves = std::move(moves);
                    }
                };
                for (int line = 0; line < config.line_search_steps; ++line) {
                    const Real scale = std::ldexp(1.0, -line);
                    const Real offset_x = block.direction_x == 0.0 ? 0.0 :
                        -std::copysign(config.step_bins * density.bin_width() * scale,
                                      block.direction_x);
                    const Real offset_y = block.direction_y == 0.0 ? 0.0 :
                        -std::copysign(config.step_bins * density.bin_height() * scale,
                                      block.direction_y);
                    consider_block(offset_x, offset_y);
                    if (config.axis_separated && offset_x != 0.0 &&
                        offset_y != 0.0) {
                        consider_block(offset_x, 0.0);
                        consider_block(0.0, offset_y);
                    }
                }
                if (config.net_block_density_direction &&
                    (block.density_direction_x != 0.0 ||
                     block.density_direction_y != 0.0)) {
                    for (int line = 0; line < config.line_search_steps; ++line) {
                        const Real scale = std::ldexp(1.0, -line);
                        const Real offset_x =
                            block.density_direction_x == 0.0 ? 0.0 :
                            -std::copysign(
                                config.step_bins * density.bin_width() * scale,
                                block.density_direction_x);
                        const Real offset_y =
                            block.density_direction_y == 0.0 ? 0.0 :
                            -std::copysign(
                                config.step_bins * density.bin_height() * scale,
                                block.density_direction_y);
                        consider_block(offset_x, offset_y);
                        if (config.axis_separated && offset_x != 0.0 &&
                            offset_y != 0.0) {
                            consider_block(offset_x, 0.0);
                            consider_block(0.0, offset_y);
                        }
                    }
                }
                if (config.breakpoint_oracle) {
                    std::vector<Real> x_offsets;
                    std::vector<Real> y_offsets;
                    const Real direction_x = block.direction_x == 0.0 ? 0.0 :
                        -std::copysign(1.0, block.direction_x);
                    const Real direction_y = block.direction_y == 0.0 ? 0.0 :
                        -std::copysign(1.0, block.direction_y);
                    for (int id : block.nodes) {
                        const Node& node = db.nodes[id];
                        for (Real x : breakpoint_positions(
                                 node.x, node.width, db.xl, db.xh,
                                 density.bin_width(), direction_x)) {
                            x_offsets.push_back(x - node.x);
                        }
                        for (Real y : breakpoint_positions(
                                 node.y, node.height, db.yl, db.yh,
                                 density.bin_height(), direction_y)) {
                            y_offsets.push_back(y - node.y);
                        }
                    }
                    auto prepare_offsets = [](std::vector<Real>& offsets) {
                        std::sort(offsets.begin(), offsets.end(),
                                  [](Real a, Real b) {
                                      return std::abs(a) < std::abs(b);
                                  });
                        offsets.erase(std::unique(
                            offsets.begin(), offsets.end(),
                            [](Real a, Real b) {
                                return std::abs(a - b) <= 1.0e-8;
                            }), offsets.end());
                        if (offsets.size() > 8) offsets.resize(8);
                    };
                    prepare_offsets(x_offsets);
                    prepare_offsets(y_offsets);
                    for (Real offset_x : x_offsets) {
                        consider_block(offset_x, 0.0);
                    }
                    for (Real offset_y : y_offsets) {
                        consider_block(0.0, offset_y);
                    }
                    const std::size_t nx = std::min<std::size_t>(4, x_offsets.size());
                    const std::size_t ny = std::min<std::size_t>(4, y_offsets.size());
                    for (std::size_t ix = 0; ix < nx; ++ix) {
                        for (std::size_t iy = 0; iy < ny; ++iy) {
                            consider_block(x_offsets[ix], y_offsets[iy]);
                        }
                    }
                }
                if (config.net_block_contraction) {
                    Real max_ratio = 1.0;
                    for (int id : block.nodes) {
                        const Node& node = db.nodes[id];
                        max_ratio = std::max(max_ratio,
                            std::abs(block.target_x - node.x) /
                                std::max<Real>(config.step_bins *
                                    density.bin_width(), 1.0e-30));
                        max_ratio = std::max(max_ratio,
                            std::abs(block.target_y - node.y) /
                                std::max<Real>(config.step_bins *
                                    density.bin_height(), 1.0e-30));
                    }
                    for (int line = 0; line < config.line_search_steps; ++line) {
                        const Real alpha = std::ldexp(1.0 / max_ratio, -line);
                        consider_contraction(alpha, 0);
                        if (config.axis_separated) {
                            consider_contraction(alpha, 1);
                            consider_contraction(alpha, 2);
                        }
                    }
                }
                if (!found) continue;
                density.commit_move(best_density_move);
                current_overflow_area += best_density_move.overflow_area_delta;
                if (!best_contraction_moves.empty()) {
                    for (const DensityNodeMove& move : best_contraction_moves) {
                        db.nodes[move.node_id].x = move.x;
                        db.nodes[move.node_id].y = move.y;
                    }
                    ++stats.net_block_contraction_moves;
                    stats.net_block_contraction_nodes +=
                        static_cast<int>(best_contraction_moves.size());
                } else {
                    for (int id : block.nodes) {
                        db.nodes[id].x += best_offset_x;
                        db.nodes[id].y += best_offset_y;
                    }
                }
                ++stats.net_block_moves;
                stats.net_block_nodes += static_cast<int>(block.nodes.size());
                ++sweep_block_moves;
            }
        }
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return std::abs(grad_x[a]) + std::abs(grad_y[a]) >
                   std::abs(grad_x[b]) + std::abs(grad_y[b]);
        });
        int sweep_moves = 0;
        const std::vector<int>& node_order = config.node_moves ? order : no_nodes;
        for (int id : node_order) {
            const bool zero_gradient =
                std::abs(grad_x[id]) <= 1.0e-18 &&
                std::abs(grad_y[id]) <= 1.0e-18;
            if (zero_gradient && !config.dead_zone_escape &&
                !config.boundary_active && !config.controlled_exact_noise) continue;
            bool bounded_degree = true;
            for (int net_id : node_nets[id]) {
                if (db.nets[net_id].pin_count >
                    static_cast<std::size_t>(config.degree_limit)) {
                    bounded_degree = false;
                    break;
                }
            }
            if (!bounded_degree) continue;

            Node& node = db.nodes[id];
            const Real direction_x = grad_x[id] +
                config.net_batch_weight * batch_x[id];
            const Real direction_y = grad_y[id] +
                config.net_batch_weight * batch_y[id];
            bool found = false;
            Real best_x = node.x;
            Real best_y = node.y;
            Real best_hpwl_delta = 0.0;
            Real best_objective_delta = 0.0;
            bool best_is_compact = false;
            DensityMove best_density_move;
            auto consider_candidate = [&](Real candidate_x, Real candidate_y,
                                          bool compact_candidate = false) {
                if (candidate_x == node.x && candidate_y == node.y) return;
                const Real hpwl_delta = affected_hpwl_delta(
                    db, id, candidate_x, candidate_y, node_nets[id]);
                DensityMove density_move = density.evaluate_move(
                    id, candidate_x, candidate_y);
                const Real objective_delta = combined_delta(
                    hpwl_delta, density_move.overflow_area_delta);
                if (objective_delta >= -1.0e-12) return;
                if (absolute_overflow_cap) {
                    if (current_overflow_area > overflow_area_cap + overflow_tolerance) {
                        // Before reaching the requested cap, allow only exact
                        // overlap-nonincreasing moves so an infeasible seed
                        // can enter the cap without introducing a surrogate.
                        if (density_move.overflow_area_delta > overflow_tolerance) {
                            return;
                        }
                    } else if (current_overflow_area + density_move.overflow_area_delta >
                               overflow_area_cap + overflow_tolerance) {
                        return;
                    }
                } else if (density_move.overflow_area_delta > overflow_tolerance) {
                    return;
                }
                if (!found || objective_delta < best_objective_delta) {
                    found = true;
                    best_x = candidate_x;
                    best_y = candidate_y;
                    best_hpwl_delta = hpwl_delta;
                    best_objective_delta = objective_delta;
                    best_is_compact = compact_candidate;
                    best_density_move = std::move(density_move);
                }
            };
            for (int line = 0; line < config.line_search_steps; ++line) {
                const Real scale = std::ldexp(1.0, -line);
                const Real step_x = config.step_bins * density.bin_width() * scale;
                const Real step_y = config.step_bins * density.bin_height() * scale;
                const Real candidate_x = std::clamp(
                    node.x - (direction_x == 0.0 ? 0.0 : std::copysign(step_x, direction_x)),
                    db.xl + 0.5 * node.width, db.xh - 0.5 * node.width);
                const Real candidate_y = std::clamp(
                    node.y - (direction_y == 0.0 ? 0.0 : std::copysign(step_y, direction_y)),
                    db.yl + 0.5 * node.height, db.yh - 0.5 * node.height);
                consider_candidate(candidate_x, candidate_y);
                if (config.axis_separated && direction_x != 0.0 && direction_y != 0.0) {
                    consider_candidate(candidate_x, node.y);
                    consider_candidate(node.x, candidate_y);
                }
            }
            if (config.breakpoint_oracle) {
                const Real breakpoint_x = direction_x == 0.0 ? 0.0 :
                    -std::copysign(1.0, direction_x);
                const Real breakpoint_y = direction_y == 0.0 ? 0.0 :
                    -std::copysign(1.0, direction_y);
                const auto x_breaks = breakpoint_positions(
                    node.x, node.width, db.xl, db.xh, density.bin_width(), breakpoint_x);
                const auto y_breaks = breakpoint_positions(
                    node.y, node.height, db.yl, db.yh, density.bin_height(), breakpoint_y);
                for (Real candidate_x : x_breaks) {
                    consider_candidate(candidate_x, node.y);
                }
                for (Real candidate_y : y_breaks) {
                    consider_candidate(node.x, candidate_y);
                }
                for (Real candidate_x : x_breaks) {
                    for (Real candidate_y : y_breaks) {
                        consider_candidate(candidate_x, candidate_y);
                    }
                }
            }
            // Dead-zone escape is an optimizer-level finite candidate oracle:
            // it is only activated when the exact first-order HPWL direction
            // is zero.  The objective itself remains exact and unsmoothed.
            if (zero_gradient && config.dead_zone_escape) {
                const auto xp = breakpoint_positions(
                    node.x, node.width, db.xl, db.xh, density.bin_width(), 1.0);
                const auto xm = breakpoint_positions(
                    node.x, node.width, db.xl, db.xh, density.bin_width(), -1.0);
                const auto yp = breakpoint_positions(
                    node.y, node.height, db.yl, db.yh, density.bin_height(), 1.0);
                const auto ym = breakpoint_positions(
                    node.y, node.height, db.yl, db.yh, density.bin_height(), -1.0);
                for (Real x : xp) consider_candidate(x, node.y);
                for (Real x : xm) consider_candidate(x, node.y);
                for (Real y : yp) consider_candidate(node.x, y);
                for (Real y : ym) consider_candidate(node.x, y);
            }
            // Boundary-strip active set: epsilon is used solely to select
            // nearby breakpoints; every proposed position is audited exactly.
            if (config.boundary_active &&
                boundary_strip_active(node.x, node.width, db.xl, db.xh,
                                      density.bin_width(), config.boundary_epsilon_bins)) {
                const auto xp = breakpoint_positions(
                    node.x, node.width, db.xl, db.xh, density.bin_width(), 1.0);
                const auto xm = breakpoint_positions(
                    node.x, node.width, db.xl, db.xh, density.bin_width(), -1.0);
                for (Real x : xp) consider_candidate(x, node.y);
                for (Real x : xm) consider_candidate(x, node.y);
            }
            if (config.boundary_active &&
                boundary_strip_active(node.y, node.height, db.yl, db.yh,
                                      density.bin_height(), config.boundary_epsilon_bins)) {
                const auto yp = breakpoint_positions(
                    node.y, node.height, db.yl, db.yh, density.bin_height(), 1.0);
                const auto ym = breakpoint_positions(
                    node.y, node.height, db.yl, db.yh, density.bin_height(), -1.0);
                for (Real y : yp) consider_candidate(node.x, y);
                for (Real y : ym) consider_candidate(node.x, y);
            }
            // Controlled exact noise: deterministic direction selection,
            // restricted to zero-gradient nodes and nearest breakpoints.
            if (zero_gradient && config.controlled_exact_noise) {
                const std::uint64_t h = mix_noise(
                    config.exact_noise_seed ^ static_cast<std::uint64_t>(id) ^
                    (static_cast<std::uint64_t>(sweep) << 32));
                const bool x_first = (h & 1ULL) != 0ULL;
                const bool positive = (h & 2ULL) != 0ULL;
                if (x_first) {
                    const auto& xs = positive
                        ? breakpoint_positions(node.x, node.width, db.xl, db.xh,
                                              density.bin_width(), 1.0)
                        : breakpoint_positions(node.x, node.width, db.xl, db.xh,
                                              density.bin_width(), -1.0);
                    for (Real x : xs) { consider_candidate(x, node.y); break; }
                } else {
                    const auto& ys = positive
                        ? breakpoint_positions(node.y, node.height, db.yl, db.yh,
                                              density.bin_height(), 1.0)
                        : breakpoint_positions(node.y, node.height, db.yl, db.yh,
                                              density.bin_height(), -1.0);
                    for (Real y : ys) { consider_candidate(node.x, y); break; }
                }
            }
            if (config.compact_directions) {
                auto consider_toward = [&](Real target_x, Real target_y) {
                    const Real dx = target_x - node.x;
                    const Real dy = target_y - node.y;
                    const Real length = std::hypot(dx, dy);
                    if (length <= 1.0e-12) return;
                    const Real max_x = config.step_bins * density.bin_width();
                    const Real max_y = config.step_bins * density.bin_height();
                    const Real normalization = std::max<Real>(
                        std::abs(dx) / std::max<Real>(max_x, 1.0e-30),
                        std::abs(dy) / std::max<Real>(max_y, 1.0e-30));
                    const Real bounded_scale = normalization > 1.0
                        ? 1.0 / normalization : 1.0;
                    for (int line = 0; line < config.line_search_steps; ++line) {
                        const Real alpha = bounded_scale * std::ldexp(1.0, -line);
                        const Real candidate_x = std::clamp(
                            node.x + alpha * dx,
                            db.xl + 0.5 * node.width,
                            db.xh - 0.5 * node.width);
                        const Real candidate_y = std::clamp(
                            node.y + alpha * dy,
                            db.yl + 0.5 * node.height,
                            db.yh - 0.5 * node.height);
                        ++stats.compact_direction_attempts;
                        consider_candidate(candidate_x, candidate_y, true);
                        if (config.axis_separated && dx != 0.0 && dy != 0.0) {
                            stats.compact_direction_attempts += 2;
                            consider_candidate(candidate_x, node.y, true);
                            consider_candidate(node.x, candidate_y, true);
                        }
                    }
                };
                consider_toward(compact_x[id], compact_y[id]);
                consider_toward(support_x, support_y);
            }
            if (!found) continue;
            density.commit_move(best_density_move);
            current_overflow_area += best_density_move.overflow_area_delta;
            node.x = best_x;
            node.y = best_y;
            ++sweep_moves;
            ++stats.moves;
            if (best_is_compact) ++stats.compact_direction_moves;
        }

        density_metrics = density.evaluate(0.0, 1.0, nullptr, nullptr);
        current_overflow_area = density_metrics.overflow * db.movable_area;
        const Real exact_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        stats.objective_evaluations += 2;
        stats.trajectory.push_back({sweep, sweep_moves + sweep_block_moves,
                                    exact_hpwl,
                                    density_metrics.overflow});
        ++stats.sweeps;
        if (sweep_moves == 0 && sweep_block_moves == 0) break;
    }

    if (stats.trajectory.empty()) {
        stats.final_hpwl = stats.initial_hpwl;
        stats.final_overflow = stats.initial_overflow;
    } else {
        stats.final_hpwl = stats.trajectory.back().hpwl;
        stats.final_overflow = stats.trajectory.back().overflow;
    }
    stats.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return stats;
}

}  // namespace ea
