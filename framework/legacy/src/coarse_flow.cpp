#include "epsilon_active/coarse_flow.hpp"

#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ea {
namespace {

struct FlowEdge {
    int source = -1;
    int destination = -1;
    Real area = 0.0;
};

struct NodeCandidate {
    Real hpwl_delta = 0.0;
    Real x = 0.0;
    Real y = 0.0;
    int id = -1;
};

struct CommodityProposal {
    std::vector<int> group;
    std::vector<DensityNodeMove> nodes;
    DensityMove density_move;
    Real hpwl_delta = 0.0;
    Real score = 0.0;
};

int coarse_bin(const Database& db, int bins_x, int bins_y, Real x, Real y) {
    const Real pitch_x = (db.xh - db.xl) / bins_x;
    const Real pitch_y = (db.yh - db.yl) / bins_y;
    const int bx = std::clamp(
        static_cast<int>(std::floor((x - db.xl) / pitch_x)), 0, bins_x - 1);
    const int by = std::clamp(
        static_cast<int>(std::floor((y - db.yl) / pitch_y)), 0, bins_y - 1);
    return by * bins_x + bx;
}

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

Real incident_hpwl_delta(const Database& db, int node_id, Real new_x, Real new_y,
                         const std::vector<int>& nets) {
    Real result = 0.0;
    for (int net_id : nets) {
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
            const Real candidate_x =
                (pin.node == node_id ? new_x : node.x) + pin.offset_x;
            const Real candidate_y =
                (pin.node == node_id ? new_y : node.y) + pin.offset_y;
            old_min_x = std::min(old_min_x, old_x);
            old_max_x = std::max(old_max_x, old_x);
            old_min_y = std::min(old_min_y, old_y);
            old_max_y = std::max(old_max_y, old_y);
            new_min_x = std::min(new_min_x, candidate_x);
            new_max_x = std::max(new_max_x, candidate_x);
            new_min_y = std::min(new_min_y, candidate_y);
            new_max_y = std::max(new_max_y, candidate_y);
        }
        result += net.weight *
            ((new_max_x - new_min_x) + (new_max_y - new_min_y) -
             (old_max_x - old_min_x) - (old_max_y - old_min_y));
    }
    return result;
}

std::vector<Real> coarse_capacity(const ExactOverlapDensity& density,
                                  int coarse_x, int coarse_y) {
    std::vector<Real> result(static_cast<std::size_t>(coarse_x) * coarse_y, 0.0);
    const Real fine_capacity = density.target_density() * density.bin_area();
    for (int by = 0; by < density.bins_y(); ++by) {
        const int cy = by * coarse_y / density.bins_y();
        for (int bx = 0; bx < density.bins_x(); ++bx) {
            const int cx = bx * coarse_x / density.bins_x();
            const int fine = by * density.bins_x() + bx;
            result[cy * coarse_x + cx] += std::max<Real>(
                fine_capacity - density.fixed_occupancy()[fine], 0.0);
        }
    }
    return result;
}

std::vector<Real> coarse_load(const Database& db, int bins_x, int bins_y,
                              std::vector<std::vector<int>>* nodes) {
    std::vector<Real> result(static_cast<std::size_t>(bins_x) * bins_y, 0.0);
    if (nodes) nodes->assign(result.size(), {});
    const Real width = (db.xh - db.xl) / bins_x;
    const Real height = (db.yh - db.yl) / bins_y;
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real left = node.x - 0.5 * node.width;
        const Real right = node.x + 0.5 * node.width;
        const Real bottom = node.y - 0.5 * node.height;
        const Real top = node.y + 0.5 * node.height;
        const int bx0 = std::clamp(static_cast<int>(std::floor(
            (left - db.xl) / width)), 0, bins_x - 1);
        const int bx1 = std::clamp(static_cast<int>(std::floor(
            (std::nextafter(right, left) - db.xl) / width)), 0, bins_x - 1);
        const int by0 = std::clamp(static_cast<int>(std::floor(
            (bottom - db.yl) / height)), 0, bins_y - 1);
        const int by1 = std::clamp(static_cast<int>(std::floor(
            (std::nextafter(top, bottom) - db.yl) / height)), 0, bins_y - 1);
        for (int by = by0; by <= by1; ++by) {
            const Real bin_bottom = db.yl + by * height;
            const Real overlap_y = std::max<Real>(
                0.0, std::min(top, bin_bottom + height) -
                      std::max(bottom, bin_bottom));
            for (int bx = bx0; bx <= bx1; ++bx) {
                const Real bin_left = db.xl + bx * width;
                const Real overlap_x = std::max<Real>(
                    0.0, std::min(right, bin_left + width) -
                      std::max(left, bin_left));
                result[by * bins_x + bx] += overlap_x * overlap_y;
            }
        }
        if (nodes) {
            const int center = coarse_bin(db, bins_x, bins_y, node.x, node.y);
            (*nodes)[center].push_back(id);
        }
    }
    return result;
}

Real group_hpwl_delta(const Database& db, const std::vector<int>& group,
                      Real dx, Real dy,
                      const std::vector<int>& incident_nets) {
    std::vector<unsigned char> moved(db.nodes.size(), 0);
    for (int id : group) moved[id] = 1;
    Real delta = 0.0;
    for (int net_id : incident_nets) {
        const Net& net = db.nets[net_id];
        if (net.pin_count < 2) continue;
        Real ominx = std::numeric_limits<Real>::infinity();
        Real omaxx = -std::numeric_limits<Real>::infinity();
        Real ominy = std::numeric_limits<Real>::infinity();
        Real omaxy = -std::numeric_limits<Real>::infinity();
        Real nminx = ominx, nmaxx = omaxx, nminy = ominy, nmaxy = omaxy;
        for (std::size_t p = net.pin_begin;
             p < net.pin_begin + net.pin_count; ++p) {
            const Pin& pin = db.pins[p];
            const Node& node = db.nodes[pin.node];
            const Real ox = node.x + pin.offset_x;
            const Real oy = node.y + pin.offset_y;
            const Real px = ox + (moved[pin.node] ? dx : 0.0);
            const Real py = oy + (moved[pin.node] ? dy : 0.0);
            ominx = std::min(ominx, ox); omaxx = std::max(omaxx, ox);
            ominy = std::min(ominy, oy); omaxy = std::max(omaxy, oy);
            nminx = std::min(nminx, px); nmaxx = std::max(nmaxx, px);
            nminy = std::min(nminy, py); nmaxy = std::max(nmaxy, py);
        }
        delta += net.weight * ((nmaxx - nminx) + (nmaxy - nminy) -
                               (omaxx - ominx) - (omaxy - ominy));
    }
    return delta;
}

Real coarse_overflow(const std::vector<Real>& load,
                     const std::vector<Real>& capacity, Real movable_area) {
    Real excess = 0.0;
    for (std::size_t i = 0; i < load.size(); ++i) {
        excess += std::max<Real>(load[i] - capacity[i], 0.0);
    }
    return excess / std::max<Real>(movable_area, 1.0e-30);
}

std::vector<FlowEdge> build_flow_plan(const std::vector<Real>& load,
                                      const std::vector<Real>& capacity,
                                      int bins_x, int bins_y, Real tolerance) {
    std::vector<Real> supply(load.size(), 0.0);
    std::vector<Real> residual(load.size(), 0.0);
    std::vector<int> sources;
    for (int bin = 0; bin < static_cast<int>(load.size()); ++bin) {
        supply[bin] = std::max<Real>(load[bin] - capacity[bin], 0.0);
        residual[bin] = std::max<Real>(capacity[bin] - load[bin], 0.0);
        if (supply[bin] > tolerance) sources.push_back(bin);
    }

    std::vector<FlowEdge> result;
    const int maximum_distance = bins_x + bins_y - 2;
    for (int distance = 1; distance <= maximum_distance; ++distance) {
        bool unfinished = false;
        for (int source : sources) {
            if (supply[source] <= tolerance) continue;
            unfinished = true;
            const int sx = source % bins_x;
            const int sy = source / bins_x;
            for (int dy = -distance; dy <= distance && supply[source] > tolerance;
                 ++dy) {
                const int dx = distance - std::abs(dy);
                const int y = sy + dy;
                if (y < 0 || y >= bins_y) continue;
                const int choices = dx == 0 ? 1 : 2;
                for (int side = 0; side < choices && supply[source] > tolerance;
                     ++side) {
                    const int x = sx + (side == 0 ? dx : -dx);
                    if (x < 0 || x >= bins_x) continue;
                    const int destination = y * bins_x + x;
                    if (residual[destination] <= tolerance) continue;
                    const Real area = std::min(supply[source], residual[destination]);
                    result.push_back({source, destination, area});
                    supply[source] -= area;
                    residual[destination] -= area;
                }
            }
        }
        if (!unfinished) break;
    }
    return result;
}

}  // namespace

CoarseFlowStats coarse_capacity_flow(
    Database& db, ExactOverlapDensity& density, const CoarseFlowConfig& config) {
    if (config.bins_x <= 0 || config.bins_y <= 0 || config.passes < 0 ||
        config.bins_x > density.bins_x() || config.bins_y > density.bins_y() ||
        config.max_source_bins < 0 || config.max_nodes_per_source < 0 ||
        config.anchor_grid_stride <= 0 || config.commodity_degree_limit < 2 ||
        config.commodity_max_nodes < 1 || config.joint_proposals_per_edge < 1 ||
        config.joint_hpwl_budget_ratio < 0.0) {
        throw std::invalid_argument("invalid coarse capacity flow configuration");
    }
    CoarseFlowStats stats;
    if (!config.enabled || config.passes == 0) return stats;

    const auto started = std::chrono::steady_clock::now();
    ExactHpwl hpwl(db);
    const DensityMetrics initial_density =
        density.evaluate(0.0, 1.0, nullptr, nullptr);
    stats.initial_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    stats.initial_overflow = initial_density.overflow;
    const std::vector<Real> capacity = coarse_capacity(
        density, config.bins_x, config.bins_y);
    const auto node_nets = build_node_nets(db);
    std::vector<Real> load = coarse_load(
        db, config.bins_x, config.bins_y, nullptr);
    stats.initial_coarse_overflow = coarse_overflow(
        load, capacity, db.movable_area);
    const Real coarse_width = (db.xh - db.xl) / config.bins_x;
    const Real coarse_height = (db.yh - db.yl) / config.bins_y;
    const Real tolerance = 1.0e-12 * coarse_width * coarse_height;

    for (int pass = 0; pass < config.passes; ++pass) {
        std::vector<std::vector<int>> nodes;
        load = coarse_load(db, config.bins_x, config.bins_y, &nodes);
        const std::vector<FlowEdge> plan = build_flow_plan(
            load, capacity, config.bins_x, config.bins_y, tolerance);
        if (plan.empty()) break;
        std::vector<FlowEdge> bounded_plan = plan;
        if (config.max_source_bins > 0) {
            std::vector<int> source_bins;
            for (const FlowEdge& edge : bounded_plan) {
                if (std::find(source_bins.begin(), source_bins.end(), edge.source) ==
                    source_bins.end()) {
                    source_bins.push_back(edge.source);
                }
            }
            std::sort(source_bins.begin(), source_bins.end(),
                      [&](int a, int b) { return load[a] - capacity[a] >
                                                 load[b] - capacity[b]; });
            if (static_cast<int>(source_bins.size()) > config.max_source_bins) {
                source_bins.resize(config.max_source_bins);
                bounded_plan.erase(std::remove_if(
                    bounded_plan.begin(), bounded_plan.end(), [&](const FlowEdge& edge) {
                        return std::find(source_bins.begin(), source_bins.end(),
                                         edge.source) == source_bins.end();
                    }), bounded_plan.end());
            }
        }
        stats.flow_edges += static_cast<int>(bounded_plan.size());
        for (const FlowEdge& edge : bounded_plan) stats.planned_area += edge.area;
        std::vector<unsigned char> moved(db.nodes.size(), 0);
        int pass_moves = 0;

        std::vector<std::vector<std::vector<int>>> commodities(nodes.size());
        if (config.connected_commodities) {
            std::vector<int> source_of(db.nodes.size(), -1);
            for (int source = 0; source < static_cast<int>(nodes.size()); ++source) {
                for (int id : nodes[source]) source_of[id] = source;
            }
            for (int source = 0; source < static_cast<int>(nodes.size()); ++source) {
                std::vector<unsigned char> assigned(db.nodes.size(), 0);
                for (int seed : nodes[source]) {
                    if (assigned[seed]) continue;
                    std::vector<int> best_group;
                    for (int net_id : node_nets[seed]) {
                        const Net& net = db.nets[net_id];
                        if (net.pin_count > static_cast<std::size_t>(
                                config.commodity_degree_limit)) continue;
                        std::vector<int> group;
                        for (std::size_t p = net.pin_begin;
                             p < net.pin_begin + net.pin_count; ++p) {
                            const int id = db.pins[p].node;
                            if (!db.nodes[id].fixed && source_of[id] == source &&
                                !assigned[id]) group.push_back(id);
                        }
                        if (group.size() >= 2 &&
                            group.size() <= static_cast<std::size_t>(
                                config.commodity_max_nodes) &&
                            group.size() > best_group.size()) {
                            best_group = std::move(group);
                        }
                    }
                    if (best_group.empty()) best_group.push_back(seed);
                    for (int id : best_group) assigned[id] = 1;
                    commodities[source].push_back(std::move(best_group));
                }
                stats.commodity_groups +=
                    static_cast<int>(commodities[source].size());
            }
        }

        if (config.connected_commodities && config.joint_assignment) {
            // Generate a small bid set per flow edge, then audit a disjoint
            // commodity batch against the exact global HPWL and overlap
            // objectives.  No proxy density or smoothed wirelength is used.
            std::vector<CommodityProposal> proposals;
            for (const FlowEdge& edge : bounded_plan) {
                const int sx = edge.source % config.bins_x;
                const int sy = edge.source / config.bins_x;
                const int dx_bin = edge.destination % config.bins_x;
                const int dy_bin = edge.destination / config.bins_x;
                const int step_x = config.exact_anchors
                    ? sx + ((dx_bin > sx) - (dx_bin < sx)) : dx_bin;
                const int step_y = config.exact_anchors
                    ? sy + ((dy_bin > sy) - (dy_bin < sy)) : dy_bin;
                const Real dx = (step_x - sx) * coarse_width;
                const Real dy = (step_y - sy) * coarse_height;
                const Real scales[] = {1.0, 0.5, 0.25};
                const Real direction_x[] = {1.0, 1.0, 0.0};
                const Real direction_y[] = {1.0, 0.0, 1.0};
                std::vector<CommodityProposal> edge_proposals;
                for (const auto& group : commodities[edge.source]) {
                    Real area = 0.0;
                    bool already_moved = false;
                    for (int id : group) {
                        area += db.nodes[id].area();
                        already_moved = already_moved || moved[id];
                    }
                    if (already_moved ||
                        (!config.joint_ignore_edge_area && area > edge.area + tolerance)) {
                        continue;
                    }
                    std::vector<int> incident_nets;
                    for (int id : group) {
                        incident_nets.insert(incident_nets.end(),
                                             node_nets[id].begin(), node_nets[id].end());
                    }
                    std::sort(incident_nets.begin(), incident_nets.end());
                    incident_nets.erase(std::unique(incident_nets.begin(),
                                                    incident_nets.end()),
                                        incident_nets.end());
                    for (int di = 0; di < 3; ++di) {
                        for (Real scale : scales) {
                            const Real sdx = dx * direction_x[di] * scale;
                            const Real sdy = dy * direction_y[di] * scale;
                            if (std::abs(sdx) < tolerance && std::abs(sdy) < tolerance) {
                                continue;
                            }
                            std::vector<DensityNodeMove> trial_nodes;
                            trial_nodes.reserve(group.size());
                            for (int id : group) {
                                const Node& node = db.nodes[id];
                                trial_nodes.push_back({id,
                                    std::clamp(node.x + sdx,
                                               db.xl + 0.5 * node.width,
                                               db.xh - 0.5 * node.width),
                                    std::clamp(node.y + sdy,
                                               db.yl + 0.5 * node.height,
                                               db.yh - 0.5 * node.height)});
                            }
                            DensityMove trial = density.evaluate_group_move(trial_nodes);
                            ++stats.commodity_scale_trials;
                            if (trial.overflow_area_delta >= -tolerance) continue;
                            const Real hd = group_hpwl_delta(
                                db, group, sdx, sdy, incident_nets);
                            if (config.commodity_hpwl_guard && hd > tolerance) continue;
                            const Real score = -trial.overflow_area_delta /
                                (std::max<Real>(hd, 0.0) +
                                 1.0e-9 * std::max<Real>(db.movable_area, 1.0));
                            edge_proposals.push_back({group, std::move(trial_nodes),
                                                      std::move(trial), hd, score});
                        }
                    }
                }
                std::stable_sort(edge_proposals.begin(), edge_proposals.end(),
                    [](const CommodityProposal& a, const CommodityProposal& b) {
                        if (a.score != b.score) return a.score > b.score;
                        return a.hpwl_delta < b.hpwl_delta;
                    });
                if (edge_proposals.size() > static_cast<std::size_t>(
                        config.joint_proposals_per_edge)) {
                    edge_proposals.resize(static_cast<std::size_t>(
                        config.joint_proposals_per_edge));
                }
                proposals.insert(proposals.end(),
                                 std::make_move_iterator(edge_proposals.begin()),
                                 std::make_move_iterator(edge_proposals.end()));
            }
            std::stable_sort(proposals.begin(), proposals.end(),
                [](const CommodityProposal& a, const CommodityProposal& b) {
                    if (a.score != b.score) return a.score > b.score;
                    return a.hpwl_delta < b.hpwl_delta;
                });
            std::vector<unsigned char> node_used(db.nodes.size(), 0);
            std::vector<CommodityProposal*> selected;
            std::vector<DensityNodeMove> batch_nodes;
            DensityMove best_batch_move;
            Real best_batch_score = -std::numeric_limits<Real>::infinity();
            Real best_batch_overflow = 0.0;
            Real best_batch_hpwl = std::numeric_limits<Real>::infinity();
            int best_count = 0;
            const Real base_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
            for (CommodityProposal& proposal : proposals) {
                if (selected.size() >= 64) break;
                bool conflict = false;
                for (int id : proposal.group) conflict = conflict || node_used[id];
                if (conflict) continue;
                selected.push_back(&proposal);
                for (const DensityNodeMove& move : proposal.nodes) {
                    node_used[move.node_id] = 1;
                    batch_nodes.push_back(move);
                }
                DensityMove batch_move = density.evaluate_group_move(batch_nodes);
                std::vector<Real> saved_x(batch_nodes.size());
                std::vector<Real> saved_y(batch_nodes.size());
                for (std::size_t mi = 0; mi < batch_nodes.size(); ++mi) {
                    const DensityNodeMove& move = batch_nodes[mi];
                    saved_x[mi] = db.nodes[move.node_id].x;
                    saved_y[mi] = db.nodes[move.node_id].y;
                    db.nodes[move.node_id].x = move.x;
                    db.nodes[move.node_id].y = move.y;
                }
                const Real trial_hpwl = hpwl.evaluate(
                    0.0, 1.0, -1, nullptr, nullptr);
                for (std::size_t mi = 0; mi < batch_nodes.size(); ++mi) {
                    const DensityNodeMove& move = batch_nodes[mi];
                    db.nodes[move.node_id].x = saved_x[mi];
                    db.nodes[move.node_id].y = saved_y[mi];
                }
                ++stats.joint_batch_trials;
                // A batch candidate is accepted only when exact overlap falls.
                // HPWL is used in the score, so cut-net cancellation is visible
                // to the prefix choice.
                if (batch_move.overflow_area_delta < -tolerance) {
                    const Real hd = trial_hpwl - base_hpwl;
                    const bool in_budget = config.joint_hpwl_budget_ratio <= 0.0 ||
                        hd <= config.joint_hpwl_budget_ratio *
                                  std::max<Real>(base_hpwl, 1.0);
                    const Real score = -batch_move.overflow_area_delta /
                        (std::max<Real>(hd, 0.0) +
                         1.0e-9 * std::max<Real>(db.movable_area, 1.0));
                    const bool better_budget = in_budget &&
                        (config.joint_hpwl_budget_ratio > 0.0 &&
                         (batch_move.overflow_area_delta < best_batch_overflow - tolerance ||
                          (std::abs(batch_move.overflow_area_delta - best_batch_overflow) <=
                               tolerance && hd < best_batch_hpwl)));
                    if (better_budget ||
                        (config.joint_hpwl_budget_ratio <= 0.0 && score > best_batch_score)) {
                        best_batch_score = score;
                        best_batch_overflow = batch_move.overflow_area_delta;
                        best_batch_hpwl = hd;
                        best_count = static_cast<int>(selected.size());
                        best_batch_move = std::move(batch_move);
                    }
                } else if (!config.joint_allow_nonmonotone) {
                    for (const DensityNodeMove& move : proposal.nodes) {
                        node_used[move.node_id] = 0;
                    }
                    batch_nodes.resize(batch_nodes.size() - proposal.nodes.size());
                    selected.pop_back();
                }
            }
            if (best_count > 0) {
                std::vector<DensityNodeMove> commit_nodes;
                for (int i = 0; i < best_count; ++i) {
                    commit_nodes.insert(commit_nodes.end(), selected[i]->nodes.begin(),
                                        selected[i]->nodes.end());
                }
                best_batch_move = density.evaluate_group_move(commit_nodes);
                density.commit_move(best_batch_move);
                for (const DensityNodeMove& move : commit_nodes) {
                    db.nodes[move.node_id].x = move.x;
                    db.nodes[move.node_id].y = move.y;
                    moved[move.node_id] = 1;
                    stats.moved_area += db.nodes[move.node_id].area();
                }
                stats.moves += best_count;
                stats.commodity_moves += best_count;
                ++stats.joint_batches;
                ++pass_moves;
            }
        }

        if (!(config.connected_commodities && config.joint_assignment)) {
        for (const FlowEdge& edge : bounded_plan) {
            const int sx = edge.source % config.bins_x;
            const int sy = edge.source / config.bins_x;
            const int dx_bin = edge.destination % config.bins_x;
            const int dy_bin = edge.destination / config.bins_x;
            // In exact-anchor mode use a one-bin step. A full source-to-sink
            // jump destroys net topology before the exact oracle can react.
            const int step_x = config.exact_anchors
                ? sx + ((dx_bin > sx) - (dx_bin < sx)) : dx_bin;
            const int step_y = config.exact_anchors
                ? sy + ((dy_bin > sy) - (dy_bin < sy)) : dy_bin;
            const Real dx = (step_x - sx) * coarse_width;
            const Real dy = (step_y - sy) * coarse_height;
            if (config.connected_commodities) {
                Real remaining = edge.area;
                auto& source_commodities = commodities[edge.source];
                std::vector<unsigned char> used(source_commodities.size(), 0);
                while (remaining > tolerance) {
                    int best_index = -1;
                    Real best_score = -std::numeric_limits<Real>::infinity();
                    Real best_overflow = 0.0;
                    Real best_hpwl = 0.0;
                    DensityMove best_move;
                    std::vector<DensityNodeMove> best_nodes;
                    for (int gi = 0; gi < static_cast<int>(source_commodities.size()); ++gi) {
                        if (used[gi]) continue;
                        const auto& group = source_commodities[gi];
                        Real area = 0.0;
                        bool already_moved = false;
                        for (int id : group) {
                            area += db.nodes[id].area();
                            already_moved = already_moved || moved[id];
                        }
                        if (already_moved || area > remaining + tolerance) continue;
                        std::vector<int> incident_nets;
                        for (int id : group) {
                            incident_nets.insert(incident_nets.end(),
                                                 node_nets[id].begin(),
                                                 node_nets[id].end());
                        }
                        std::sort(incident_nets.begin(), incident_nets.end());
                        incident_nets.erase(std::unique(incident_nets.begin(),
                                                        incident_nets.end()),
                                            incident_nets.end());
                        const Real scales[] = {1.0, 0.5, 0.25};
                        const int scale_count = config.commodity_multiscale ? 3 : 1;
                        const Real direction_x[] = {1.0, 1.0, 0.0};
                        const Real direction_y[] = {1.0, 0.0, 1.0};
                        const int direction_count = config.commodity_axis_split ? 3 : 1;
                        for (int di = 0; di < direction_count; ++di) {
                            if (direction_x[di] == 0.0 && direction_y[di] == 0.0) continue;
                            for (int si = 0; si < scale_count; ++si) {
                                const Real scale = scales[si];
                                const Real sdx = dx * direction_x[di] * scale;
                                const Real sdy = dy * direction_y[di] * scale;
                            std::vector<DensityNodeMove> trial_nodes;
                            trial_nodes.reserve(group.size());
                            for (int id : group) {
                                const Node& node = db.nodes[id];
                                trial_nodes.push_back({id,
                                    std::clamp(node.x + sdx,
                                               db.xl + 0.5 * node.width,
                                               db.xh - 0.5 * node.width),
                                    std::clamp(node.y + sdy,
                                               db.yl + 0.5 * node.height,
                                               db.yh - 0.5 * node.height)});
                            }
                            DensityMove trial = density.evaluate_group_move(trial_nodes);
                                ++stats.commodity_scale_trials;
                                if (trial.overflow_area_delta >= -tolerance) continue;
                                const Real hd = group_hpwl_delta(
                                    db, group, sdx, sdy, incident_nets);
                                if (config.commodity_hpwl_guard && hd > tolerance) {
                                    continue;
                                }
                                const Real score = -trial.overflow_area_delta /
                                    (std::max<Real>(hd, 0.0) +
                                     1.0e-9 * std::max<Real>(db.movable_area, 1.0));
                                if (best_index < 0 || score > best_score ||
                                    (score == best_score && hd < best_hpwl)) {
                                    best_index = gi;
                                    best_score = score;
                                    best_overflow = trial.overflow_area_delta;
                                    best_hpwl = hd;
                                    best_move = std::move(trial);
                                    best_nodes = std::move(trial_nodes);
                                }
                            }
                        }
                    }
                    if (best_index < 0) break;
                    density.commit_move(best_move);
                    used[best_index] = 1;
                    Real area = 0.0;
                    for (const DensityNodeMove& move : best_nodes) {
                        db.nodes[move.node_id].x = move.x;
                        db.nodes[move.node_id].y = move.y;
                        moved[move.node_id] = 1;
                        area += db.nodes[move.node_id].area();
                    }
                    (void)best_overflow;
                    remaining -= area;
                    stats.moved_area += area;
                    ++stats.moves;
                    ++stats.commodity_moves;
                    ++pass_moves;
                }
                continue;
            }
            std::vector<NodeCandidate> candidates;
            candidates.reserve(nodes[edge.source].size());
            for (int id : nodes[edge.source]) {
                if (moved[id]) continue;
                const Node& node = db.nodes[id];
                const Real x = std::clamp(
                    node.x + dx, db.xl + 0.5 * node.width,
                    db.xh - 0.5 * node.width);
                const Real y = std::clamp(
                    node.y + dy, db.yl + 0.5 * node.height,
                    db.yh - 0.5 * node.height);
                candidates.push_back({0.0, x, y, id});
            }
            if (config.max_nodes_per_source > 0 &&
                static_cast<int>(candidates.size()) > config.max_nodes_per_source) {
                std::stable_sort(candidates.begin(), candidates.end(),
                                 [](const NodeCandidate& a, const NodeCandidate& b) {
                                     return a.id < b.id;
                                 });
                candidates.resize(config.max_nodes_per_source);
            }
            #pragma omp parallel for schedule(dynamic, 32) if(candidates.size() > 64)
            for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
                NodeCandidate& candidate = candidates[i];
                candidate.hpwl_delta = incident_hpwl_delta(
                    db, candidate.id, candidate.x, candidate.y,
                    node_nets[candidate.id]);
            }
            std::sort(candidates.begin(), candidates.end(),
                      [](const NodeCandidate& a, const NodeCandidate& b) {
                if (a.hpwl_delta != b.hpwl_delta) {
                    return a.hpwl_delta < b.hpwl_delta;
                }
                return a.id < b.id;
            });
            Real remaining = edge.area;
            for (const NodeCandidate& candidate : candidates) {
                if (remaining <= tolerance) break;
                Node& node = db.nodes[candidate.id];
                if (config.exact_anchors) {
                    const int fine_x0 = step_x * density.bins_x() / config.bins_x;
                    const int fine_x1 = (step_x + 1) * density.bins_x() /
                        config.bins_x;
                    const int fine_y0 = step_y * density.bins_y() / config.bins_y;
                    const int fine_y1 = (step_y + 1) * density.bins_y() /
                        config.bins_y;
                    bool found = false;
                    Real best_overflow = 0.0;
                    Real best_hpwl = std::numeric_limits<Real>::infinity();
                    int best_bin = -1;
                    Real best_x = node.x;
                    Real best_y = node.y;
                    DensityMove best_move;
                    for (int fine_y = fine_y0; fine_y < fine_y1;
                         fine_y += config.anchor_grid_stride) {
                        for (int fine_x = fine_x0; fine_x < fine_x1;
                             fine_x += config.anchor_grid_stride) {
                            const Real x = std::clamp(
                                db.xl + (fine_x + 0.5) * density.bin_width(),
                                db.xl + 0.5 * node.width,
                                db.xh - 0.5 * node.width);
                            const Real y = std::clamp(
                                db.yl + (fine_y + 0.5) * density.bin_height(),
                                db.yl + 0.5 * node.height,
                                db.yh - 0.5 * node.height);
                            DensityMove move = density.evaluate_move(
                                candidate.id, x, y);
                            ++stats.exact_anchor_trials;
                            if (move.overflow_area_delta >= -tolerance) continue;
                            const Real hpwl_delta = incident_hpwl_delta(
                                db, candidate.id, x, y,
                                node_nets[candidate.id]);
                            const int fine_bin = fine_y * density.bins_x() + fine_x;
                            if (!found ||
                                move.overflow_area_delta < best_overflow - tolerance ||
                                (std::abs(move.overflow_area_delta - best_overflow) <=
                                     tolerance &&
                                 (hpwl_delta < best_hpwl - 1.0e-12 ||
                                  (std::abs(hpwl_delta - best_hpwl) <= 1.0e-12 &&
                                   fine_bin < best_bin)))) {
                                found = true;
                                best_overflow = move.overflow_area_delta;
                                best_hpwl = hpwl_delta;
                                best_bin = fine_bin;
                                best_x = x;
                                best_y = y;
                                best_move = std::move(move);
                            }
                        }
                    }
                    if (!found) {
                        ++stats.exact_anchor_rejections;
                        continue;
                    }
                    node.x = best_x;
                    node.y = best_y;
                    density.commit_move(best_move);
                } else {
                    node.x = candidate.x;
                    node.y = candidate.y;
                }
                moved[candidate.id] = 1;
                remaining -= node.area();
                stats.moved_area += node.area();
                ++stats.moves;
                ++pass_moves;
            }
        }
        }
        ++stats.passes;
        if (pass_moves == 0) break;
    }

    load = coarse_load(db, config.bins_x, config.bins_y, nullptr);
    stats.final_coarse_overflow = coarse_overflow(load, capacity, db.movable_area);
    const DensityMetrics final_density =
        density.evaluate(0.0, 1.0, nullptr, nullptr);
    stats.final_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    stats.final_overflow = final_density.overflow;
    stats.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return stats;
}

}  // namespace ea
