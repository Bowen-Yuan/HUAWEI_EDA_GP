#include "epsilon_active/density_coordinate.hpp"

#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace ea {
namespace {

Real affected_hpwl_delta(const Database& db, int node_id, Real nx, Real ny,
                         const std::vector<int>& incident) {
    Real delta = 0.0;
    for (int net_id : incident) {
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
            const Real px = (pin.node == node_id ? nx : node.x) + pin.offset_x;
            const Real py = (pin.node == node_id ? ny : node.y) + pin.offset_y;
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

Real affected_hpwl_delta_group(
    const Database& db, const std::vector<int>& nodes, Real dx, Real dy,
    const std::vector<int>& incident, std::vector<int>& marks, int mark) {
    for (int id : nodes) marks[id] = mark;
    Real delta = 0.0;
    for (int net_id : incident) {
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
            const bool moved = marks[pin.node] == mark;
            const Real px = ox + (moved ? dx : 0.0);
            const Real py = oy + (moved ? dy : 0.0);
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

// Exact HPWL delta for an arbitrary per-node move set.  This is intentionally
// a piecewise max/min oracle; no smooth surrogate is introduced.
Real affected_hpwl_delta_moves(
    const Database& db, const std::vector<DensityNodeMove>& moves,
    const std::vector<int>& incident) {
    Real delta = 0.0;
    for (int net_id : incident) {
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
            Real px = ox;
            Real py = oy;
            for (const DensityNodeMove& move : moves) {
                if (move.node_id == pin.node) {
                    px = move.x + pin.offset_x;
                    py = move.y + pin.offset_y;
                    break;
                }
            }
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

std::vector<Real> nearest_breakpoints(Real center, Real extent, Real lower,
                                       Real upper, Real pitch, Real direction) {
    std::vector<Real> result;
    if (direction == 0.0 || pitch <= 0.0) return result;
    const Real edge = direction > 0.0 ? center - 0.5 * extent
                                      : center + 0.5 * extent;
    const Real shifted = (edge - lower) / pitch;
    const Real first = direction > 0.0 ? std::floor(shifted) + 1.0
                                       : std::ceil(shifted) - 1.0;
    for (int k = 0; k < 4; ++k) {
        const Real boundary = lower +
            (first + direction * static_cast<Real>(k)) * pitch;
        const Real distance = direction * (boundary - edge);
        if (distance <= 1.0e-8) continue;
        const Real candidate = std::clamp(
            center + direction * distance,
            lower + 0.5 * extent, upper - 0.5 * extent);
        if (std::abs(candidate - center) > 1.0e-8) result.push_back(candidate);
    }
    return result;
}

}  // namespace

DensityCoordinateStats coordinate_descent_overlap(
    Database& db, ExactOverlapDensity& density,
    const DensityCoordinateConfig& config) {
    if (config.sweeps < 0 || config.line_search_steps <= 0 ||
        config.step_bins <= 0.0 || config.hpwl_budget_fraction < 0.0 ||
        config.overflow_band_fraction < 0.0 ||
        config.degree_limit < 2 || config.assignment_bins <= 0 ||
        config.hierarchy_levels <= 0 || config.hierarchy_levels > 16 ||
        config.cluster_max_nodes < 2 ||
        config.cluster_candidate_bins <= 0 ||
        config.cluster_max_clusters < 0 ||
        config.adjacent_bin_radius < 1 ||
        config.block_degree_limit < 2 ||
        config.block_max_nodes < 2 || config.block_max_blocks < 0) {
        throw std::invalid_argument("invalid overlap coordinate configuration");
    }
    DensityCoordinateStats stats;
    if (config.sweeps == 0) return stats;
    const auto started = std::chrono::steady_clock::now();
    ExactHpwl hpwl(db);
    DensityMetrics current = density.evaluate(0.0, 1.0, nullptr, nullptr);
    Real current_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    const Real initial_hpwl = current_hpwl;
    const Real budget = initial_hpwl * (1.0 + config.hpwl_budget_fraction);
    stats.initial_hpwl = current_hpwl;
    stats.initial_overflow = current.overflow;

    // Acceptance must include every incident net. Degree limits are only
    // candidate-generation controls, never objective approximations.
    std::vector<std::vector<int>> incident(db.nodes.size());
    std::vector<std::vector<int>> connectivity(db.nodes.size());
    for (int net_id = 0; net_id < static_cast<int>(db.nets.size()); ++net_id) {
        const Net& net = db.nets[net_id];
        for (std::size_t p = net.pin_begin;
             p < net.pin_begin + net.pin_count; ++p) {
            incident[db.pins[p].node].push_back(net_id);
            if (net.pin_count <=
                static_cast<std::size_t>(config.block_degree_limit)) {
                connectivity[db.pins[p].node].push_back(net_id);
            }
        }
    }

    const Real tolerance = 1.0e-12 * std::max<Real>(db.movable_area, 1.0);
    for (int sweep = 0; sweep < config.sweeps; ++sweep) {
        int sweep_moves = 0;
        std::vector<unsigned char> block_touched(db.nodes.size(), 0);
        std::vector<int> marks(db.nodes.size(), 0);
        int mark = 0;
        const int coarse_x = std::min(config.assignment_bins, density.bins_x());
        const int coarse_y = std::min(config.assignment_bins, density.bins_y());
        std::vector<Real> coarse_residual;
        if (config.net_blocks && config.capacity_assignment) {
            coarse_residual.assign(static_cast<std::size_t>(coarse_x) * coarse_y,
                                   0.0);
            const auto& occupancy = density.occupancy();
            const Real fine_capacity = density.target_density() * density.bin_area();
            for (int by = 0; by < density.bins_y(); ++by) {
                const int cy = by * coarse_y / density.bins_y();
                for (int bx = 0; bx < density.bins_x(); ++bx) {
                    const int cx = bx * coarse_x / density.bins_x();
                    const int fine = by * density.bins_x() + bx;
                    const int coarse = cy * coarse_x + cx;
                    // Total occupancy includes movable rectangles and fixed
                    // macros. Fixed area is therefore unavailable capacity.
                    coarse_residual[coarse] +=
                        fine_capacity - occupancy[fine];
                }
            }
        }

        if (config.hierarchical_clusters) {
            struct ClusterProposal {
                std::vector<int> nodes;
                std::vector<DensityNodeMove> moves;
                int destination = -1;
                Real dx = 0.0;
                Real dy = 0.0;
                Real overflow_delta = 0.0;
                Real hpwl_delta = 0.0;
                Real score = 0.0;
            };

            const int finest_x = std::min(
                config.assignment_bins, density.bins_x());
            const int finest_y = std::min(
                config.assignment_bins, density.bins_y());
            int previous_x = -1;
            int previous_y = -1;
            for (int level = 0; level < config.hierarchy_levels; ++level) {
                const int divisor = 1 <<
                    (config.hierarchy_levels - level - 1);
                const int level_x = std::max(1, finest_x / divisor);
                const int level_y = std::max(1, finest_y / divisor);
                if (level_x == previous_x && level_y == previous_y) continue;
                previous_x = level_x;
                previous_y = level_y;

                std::vector<Real> residual(
                    static_cast<std::size_t>(level_x) * level_y, 0.0);
                const auto& occupancy = density.occupancy();
                const Real fine_capacity =
                    density.target_density() * density.bin_area();
                for (int fy = 0; fy < density.bins_y(); ++fy) {
                    const int cy = fy * level_y / density.bins_y();
                    for (int fx = 0; fx < density.bins_x(); ++fx) {
                        const int cx = fx * level_x / density.bins_x();
                        const int fine = fy * density.bins_x() + fx;
                        residual[cy * level_x + cx] +=
                            fine_capacity - occupancy[fine];
                    }
                }

                const auto node_bin = [&](int id) {
                    const Node& node = db.nodes[id];
                    const int bx = std::clamp(
                        static_cast<int>((node.x - db.xl) * level_x /
                                         (db.xh - db.xl)),
                        0, level_x - 1);
                    const int by = std::clamp(
                        static_cast<int>((node.y - db.yl) * level_y /
                                         (db.yh - db.yl)),
                        0, level_y - 1);
                    return by * level_x + bx;
                };

                std::vector<std::vector<int>> source_nodes(residual.size());
                for (int id : db.movable_ids) {
                    source_nodes[node_bin(id)].push_back(id);
                }
                std::vector<int> source_order(residual.size());
                std::iota(source_order.begin(), source_order.end(), 0);
                std::sort(source_order.begin(), source_order.end(),
                          [&](int lhs, int rhs) {
                              return residual[lhs] < residual[rhs];
                          });

                std::vector<unsigned char> reserved(db.nodes.size(), 0);
                std::vector<int> queued(db.nodes.size(), 0);
                std::vector<int> included(db.nodes.size(), 0);
                std::vector<Real> frontier_score(db.nodes.size(), 0.0);
                int queue_token = 0;
                int cluster_attempts = 0;
                std::vector<ClusterProposal> proposals;
                proposals.reserve(config.cluster_max_clusters);

                for (int source : source_order) {
                    if (cluster_attempts >=
                        config.cluster_max_clusters) break;
                    if (residual[source] >= -tolerance) break;
                    auto& seeds = source_nodes[source];
                    std::sort(seeds.begin(), seeds.end(), [&](int lhs, int rhs) {
                        return db.nodes[lhs].area() > db.nodes[rhs].area();
                    });
                    Real largest_destination = 0.0;
                    for (Real value : residual) {
                        largest_destination = std::max(largest_destination, value);
                    }
                    if (largest_destination <= tolerance) break;
                    const Real target_area = std::min(
                        -residual[source], largest_destination);

                    for (int seed : seeds) {
                        if (cluster_attempts >=
                            config.cluster_max_clusters) break;
                        if (config.cooperative_batch &&
                            residual[source] >= -tolerance) break;
                        if (reserved[seed]) continue;
                        ++cluster_attempts;
                        ++queue_token;
                        std::vector<int> cluster;
                        Real cluster_area = 0.0;
                        if (config.heavy_edge_growth) {
                            std::priority_queue<std::pair<Real, int>> frontier;
                            queued[seed] = queue_token;
                            frontier_score[seed] =
                                std::numeric_limits<Real>::infinity();
                            frontier.emplace(frontier_score[seed], seed);
                            while (!frontier.empty() &&
                                   static_cast<int>(cluster.size()) <
                                       config.cluster_max_nodes) {
                                const auto [priority, id] = frontier.top();
                                frontier.pop();
                                if (reserved[id] ||
                                    included[id] == queue_token ||
                                    node_bin(id) != source) continue;
                                if (priority + 1.0e-15 < frontier_score[id])
                                    continue;
                                included[id] = queue_token;
                                cluster.push_back(id);
                                cluster_area += db.nodes[id].area();
                                if (cluster.size() >= 2 &&
                                    cluster_area >= target_area) break;
                                for (int net_id : connectivity[id]) {
                                    const Net& net = db.nets[net_id];
                                    const Real edge_weight = net.weight /
                                        std::max<Real>(
                                            static_cast<Real>(net.pin_count) - 1.0,
                                            1.0);
                                    for (std::size_t p = net.pin_begin;
                                         p < net.pin_begin + net.pin_count; ++p) {
                                        const int neighbor = db.pins[p].node;
                                        if (db.nodes[neighbor].fixed ||
                                            reserved[neighbor] ||
                                            included[neighbor] == queue_token ||
                                            node_bin(neighbor) != source) continue;
                                        if (queued[neighbor] != queue_token) {
                                            queued[neighbor] = queue_token;
                                            frontier_score[neighbor] = 0.0;
                                        }
                                        frontier_score[neighbor] += edge_weight;
                                        frontier.emplace(
                                            frontier_score[neighbor], neighbor);
                                    }
                                }
                            }
                        } else {
                            std::deque<int> frontier;
                            frontier.push_back(seed);
                            queued[seed] = queue_token;
                            while (!frontier.empty() &&
                                   static_cast<int>(cluster.size()) <
                                       config.cluster_max_nodes) {
                                const int id = frontier.front();
                                frontier.pop_front();
                                if (reserved[id] || node_bin(id) != source)
                                    continue;
                                cluster.push_back(id);
                                cluster_area += db.nodes[id].area();
                                if (cluster.size() >= 2 &&
                                    cluster_area >= target_area) break;
                                for (int net_id : connectivity[id]) {
                                    const Net& net = db.nets[net_id];
                                    for (std::size_t p = net.pin_begin;
                                         p < net.pin_begin + net.pin_count; ++p) {
                                        const int neighbor = db.pins[p].node;
                                        if (db.nodes[neighbor].fixed ||
                                            reserved[neighbor] ||
                                            queued[neighbor] == queue_token ||
                                            node_bin(neighbor) != source) continue;
                                        queued[neighbor] = queue_token;
                                        frontier.push_back(neighbor);
                                    }
                                }
                            }
                        }
                        if (cluster.size() < 2) continue;

                        Real center_x = 0.0;
                        Real center_y = 0.0;
                        Real dx_low = -std::numeric_limits<Real>::infinity();
                        Real dx_high = std::numeric_limits<Real>::infinity();
                        Real dy_low = -std::numeric_limits<Real>::infinity();
                        Real dy_high = std::numeric_limits<Real>::infinity();
                        std::vector<int> incident_nets;
                        for (int id : cluster) {
                            const Node& node = db.nodes[id];
                            center_x += node.x;
                            center_y += node.y;
                            dx_low = std::max(
                                dx_low, db.xl + 0.5 * node.width - node.x);
                            dx_high = std::min(
                                dx_high, db.xh - 0.5 * node.width - node.x);
                            dy_low = std::max(
                                dy_low, db.yl + 0.5 * node.height - node.y);
                            dy_high = std::min(
                                dy_high, db.yh - 0.5 * node.height - node.y);
                            incident_nets.insert(incident_nets.end(),
                                                 incident[id].begin(),
                                                 incident[id].end());
                        }
                        center_x /= static_cast<Real>(cluster.size());
                        center_y /= static_cast<Real>(cluster.size());
                        std::sort(incident_nets.begin(), incident_nets.end());
                        incident_nets.erase(
                            std::unique(incident_nets.begin(),
                                        incident_nets.end()),
                            incident_nets.end());

                        std::vector<std::tuple<Real, int>> destinations;
                        destinations.reserve(residual.size());
                        for (int destination = 0;
                             destination < static_cast<int>(residual.size());
                             ++destination) {
                            if (destination == source ||
                                residual[destination] <= tolerance) continue;
                            const int bx = destination % level_x;
                            const int by = destination / level_x;
                            if (config.adjacent_bin_only) {
                                const int sx = source % level_x;
                                const int sy = source / level_x;
                                const int manhattan = std::abs(bx - sx) +
                                    std::abs(by - sy);
                                if (manhattan > config.adjacent_bin_radius) continue;
                            }
                            const Real tx = db.xl + (bx + 0.5) *
                                (db.xh - db.xl) / level_x;
                            const Real ty = db.yl + (by + 0.5) *
                                (db.yh - db.yl) / level_y;
                            const Real shortfall = std::max<Real>(
                                cluster_area - residual[destination], 0.0) /
                                std::max<Real>(cluster_area, 1.0e-30);
                            const Real distance = std::hypot(
                                (tx - center_x) / (db.xh - db.xl),
                                (ty - center_y) / (db.yh - db.yl));
                            destinations.emplace_back(
                                distance + 2.0 * shortfall, destination);
                        }
                        std::sort(destinations.begin(), destinations.end());
                        const int candidate_count = std::min<int>(
                            config.cluster_candidate_bins,
                            static_cast<int>(destinations.size()));
                        bool found = false;
                        ClusterProposal best;
                        for (int candidate = 0; candidate < candidate_count;
                             ++candidate) {
                            const int destination =
                                std::get<1>(destinations[candidate]);
                            const int bx = destination % level_x;
                            const int by = destination / level_x;
                            const Real tx = db.xl + (bx + 0.5) *
                                (db.xh - db.xl) / level_x;
                            const Real ty = db.yl + (by + 0.5) *
                                (db.yh - db.yl) / level_y;
                            const Real full_dx = std::clamp(
                                tx - center_x, dx_low, dx_high);
                            const Real full_dy = std::clamp(
                                ty - center_y, dy_low, dy_high);
                            for (int line = 0;
                                 line < config.line_search_steps; ++line) {
                                const Real scale = std::ldexp(1.0, -line);
                                const Real dx = full_dx * scale;
                                const Real dy = full_dy * scale;
                                if (std::abs(dx) <= 1.0e-10 &&
                                    std::abs(dy) <= 1.0e-10) continue;
                                ++stats.candidates;
                                std::vector<DensityNodeMove> moves;
                                moves.reserve(cluster.size());
                                for (int id : cluster) {
                                    Real alpha = 1.0;
                                    if (config.elastic_clusters) {
                                        const Real fractions[] = {1.0, 0.75,
                                                                  0.5, 0.25,
                                                                  0.0};
                                        Real best_single =
                                            std::numeric_limits<Real>::infinity();
                                        for (Real fraction : fractions) {
                                            const Node& node = db.nodes[id];
                                            const Real nx = std::clamp(
                                                node.x + fraction * dx,
                                                db.xl + 0.5 * node.width,
                                                db.xh - 0.5 * node.width);
                                            const Real ny = std::clamp(
                                                node.y + fraction * dy,
                                                db.yl + 0.5 * node.height,
                                                db.yh - 0.5 * node.height);
                                            const Real single = affected_hpwl_delta(
                                                db, id, nx, ny, incident[id]);
                                            if (single < best_single) {
                                                best_single = single;
                                                alpha = fraction;
                                            }
                                        }
                                    }
                                    const Node& node = db.nodes[id];
                                    moves.push_back({id,
                                        std::clamp(node.x + alpha * dx,
                                                   db.xl + 0.5 * node.width,
                                                   db.xh - 0.5 * node.width),
                                        std::clamp(node.y + alpha * dy,
                                                   db.yl + 0.5 * node.height,
                                                   db.yh - 0.5 * node.height)});
                                }
                                const DensityMove trial =
                                    density.evaluate_group_move(moves);
                                const Real overflow_band =
                                    config.overflow_band_fraction *
                                    std::max<Real>(db.movable_area, 1.0);
                                if (trial.overflow_area_delta >= overflow_band)
                                    continue;
                                ++mark;
                                const Real hd = config.elastic_clusters
                                    ? affected_hpwl_delta_moves(
                                        db, moves, incident_nets)
                                    : affected_hpwl_delta_group(
                                        db, cluster, dx, dy, incident_nets,
                                        marks, mark);
                                if (!config.cooperative_batch &&
                                    current_hpwl + hd > budget + 1.0e-6)
                                    continue;
                                const Real hpwl_cost = std::max<Real>(hd, 0.0) +
                                    1.0e-9 * std::max<Real>(initial_hpwl, 1.0);
                                const Real score =
                                    -trial.overflow_area_delta / hpwl_cost;
                                if (!found || score > best.score ||
                                    (score == best.score &&
                                     trial.overflow_area_delta <
                                         best.overflow_delta)) {
                                    found = true;
                                    best.nodes = cluster;
                                    if (config.elastic_clusters) {
                                        best.moves = moves;
                                    }
                                    best.destination = destination;
                                    best.dx = dx;
                                    best.dy = dy;
                                    best.overflow_delta =
                                        trial.overflow_area_delta;
                                    best.hpwl_delta = hd;
                                    best.score = score;
                                }
                            }
                        }
                        if (!found) continue;
                        for (int id : cluster) reserved[id] = 1;
                        if (config.cooperative_batch) {
                            residual[best.destination] -= cluster_area;
                            residual[source] += cluster_area;
                        }
                        proposals.push_back(std::move(best));
                        ++stats.cluster_proposals;
                    }
                }

                std::sort(proposals.begin(), proposals.end(),
                          [](const ClusterProposal& lhs,
                             const ClusterProposal& rhs) {
                              if (lhs.score != rhs.score)
                                  return lhs.score > rhs.score;
                              return lhs.overflow_delta < rhs.overflow_delta;
                          });
                if (config.cooperative_batch && !proposals.empty()) {
                    const DensityMetrics baseline_density =
                        density.evaluate(0.0, 1.0, nullptr, nullptr);
                    std::size_t prefix = proposals.size();
                    while (prefix > 0) {
                        ++stats.cooperative_batch_trials;
                        for (std::size_t i = 0; i < prefix; ++i) {
                            const ClusterProposal& proposal = proposals[i];
                            for (int id : proposal.nodes) {
                                db.nodes[id].x += proposal.dx;
                                db.nodes[id].y += proposal.dy;
                            }
                        }
                        const DensityMetrics candidate_density =
                            density.evaluate(0.0, 1.0, nullptr, nullptr);
                        const Real candidate_hpwl =
                            hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
                        if (candidate_density.overflow <=
                                baseline_density.overflow +
                                    config.overflow_band_fraction &&
                            (config.overflow_band_fraction > 0.0 ||
                             candidate_density.overflow < baseline_density.overflow -
                                tolerance / std::max<Real>(db.movable_area, 1.0)) &&
                            candidate_hpwl <= budget + 1.0e-6) {
                            current_hpwl = candidate_hpwl;
                            sweep_moves += static_cast<int>(prefix);
                            stats.moves += static_cast<int>(prefix);
                            stats.cluster_moves += static_cast<int>(prefix);
                            ++stats.cooperative_batches;
                            for (std::size_t i = 0; i < prefix; ++i) {
                                for (int id : proposals[i].nodes) {
                                    block_touched[id] = 1;
                                }
                            }
                            break;
                        }
                        for (std::size_t i = 0; i < prefix; ++i) {
                            const ClusterProposal& proposal = proposals[i];
                            for (int id : proposal.nodes) {
                                db.nodes[id].x -= proposal.dx;
                                db.nodes[id].y -= proposal.dy;
                            }
                        }
                        density.evaluate(0.0, 1.0, nullptr, nullptr);
                        prefix /= 2;
                    }
                } else for (const ClusterProposal& proposal : proposals) {
                    if (residual[proposal.destination] <= tolerance) continue;
                    std::vector<int> incident_nets;
                    std::vector<DensityNodeMove> moves;
                    if (!proposal.moves.empty()) {
                        moves = proposal.moves;
                    } else {
                        moves.reserve(proposal.nodes.size());
                        for (int id : proposal.nodes) {
                            moves.push_back({id,
                                db.nodes[id].x + proposal.dx,
                                db.nodes[id].y + proposal.dy});
                        }
                    }
                    for (int id : proposal.nodes) {
                        incident_nets.insert(incident_nets.end(),
                                             incident[id].begin(),
                                             incident[id].end());
                    }
                    std::sort(incident_nets.begin(), incident_nets.end());
                    incident_nets.erase(
                        std::unique(incident_nets.begin(),
                                    incident_nets.end()),
                        incident_nets.end());
                    const DensityMove trial = density.evaluate_group_move(moves);
                    if (trial.overflow_area_delta >= -tolerance) continue;
                    ++mark;
                    const Real hd = proposal.moves.empty()
                        ? affected_hpwl_delta_group(
                            db, proposal.nodes, proposal.dx, proposal.dy,
                            incident_nets, marks, mark)
                        : affected_hpwl_delta_moves(db, moves, incident_nets);
                    if (current_hpwl + hd > budget + 1.0e-6) continue;
                    density.commit_move(trial);
                    for (const auto& [fine, delta] : trial.area_changes) {
                        const int fx = fine % density.bins_x();
                        const int fy = fine / density.bins_x();
                        const int cx = fx * level_x / density.bins_x();
                        const int cy = fy * level_y / density.bins_y();
                        residual[cy * level_x + cx] -= delta;
                    }
                    for (const DensityNodeMove& move : moves) {
                        db.nodes[move.node_id].x = move.x;
                        db.nodes[move.node_id].y = move.y;
                        block_touched[move.node_id] = 1;
                    }
                    current_hpwl += hd;
                    ++sweep_moves;
                    ++stats.moves;
                    ++stats.cluster_moves;
                }
                current_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
            }
        }

        // A block pass moves all movable pins of a small net rigidly.  The
        // candidate is still accepted by the exact overlap oracle and the
        // same cumulative HPWL guard as node moves.
        if (config.net_blocks && !config.hierarchical_clusters) {
            int block_count = 0;
            for (const Net& net : db.nets) {
                if (block_count >= config.block_max_blocks || net.pin_count < 2 ||
                    static_cast<int>(net.pin_count) > config.block_degree_limit) {
                    continue;
                }
                std::vector<int> block;
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const int id = db.pins[p].node;
                    if (!db.nodes[id].fixed) block.push_back(id);
                }
                std::sort(block.begin(), block.end());
                block.erase(std::unique(block.begin(), block.end()), block.end());
                if (block.size() < 2 || block.size() >
                    static_cast<std::size_t>(config.block_max_nodes)) continue;
                bool touched = false;
                for (int id : block) touched = touched || block_touched[id];
                if (touched) continue;
                ++block_count;
                std::vector<int> incident_nets;
                for (int id : block) {
                    incident_nets.insert(incident_nets.end(),
                                         incident[id].begin(), incident[id].end());
                }
                std::sort(incident_nets.begin(), incident_nets.end());
                incident_nets.erase(std::unique(incident_nets.begin(),
                                                incident_nets.end()),
                                    incident_nets.end());
                bool found = false;
                Real best_dx = 0.0, best_dy = 0.0;
                Real best_overflow_delta = 0.0;
                Real best_hpwl_delta = std::numeric_limits<Real>::infinity();
                auto consider_block = [&](Real dx, Real dy) {
                    if (std::abs(dx) <= 1.0e-10 && std::abs(dy) <= 1.0e-10) return;
                    ++stats.candidates;
                    std::vector<DensityNodeMove> moves;
                    moves.reserve(block.size());
                    for (int id : block) {
                        const Node& node = db.nodes[id];
                        moves.push_back({id,
                            std::clamp(node.x + dx,
                                       db.xl + 0.5 * node.width,
                                       db.xh - 0.5 * node.width),
                            std::clamp(node.y + dy,
                                       db.yl + 0.5 * node.height,
                                       db.yh - 0.5 * node.height)});
                    }
                    DensityMove trial = density.evaluate_group_move(moves);
                    if (trial.overflow_area_delta >= -tolerance) return;
                    ++mark;
                    const Real hd = affected_hpwl_delta_group(
                        db, block, dx, dy, incident_nets, marks, mark);
                    if (current_hpwl + hd > budget + 1.0e-6) return;
                    if (!found || trial.overflow_area_delta <
                                      best_overflow_delta - tolerance ||
                        (std::abs(trial.overflow_area_delta -
                                  best_overflow_delta) <= tolerance &&
                         hd < best_hpwl_delta)) {
                        found = true;
                        best_dx = dx;
                        best_dy = dy;
                        best_overflow_delta = trial.overflow_area_delta;
                        best_hpwl_delta = hd;
                    }
                };
                for (int line = 0; line < config.line_search_steps; ++line) {
                    const Real scale = std::ldexp(1.0, -line);
                    const Real dx = config.step_bins * density.bin_width() * scale;
                    const Real dy = config.step_bins * density.bin_height() * scale;
                    consider_block(-dx, 0.0); consider_block(dx, 0.0);
                    consider_block(0.0, -dy); consider_block(0.0, dy);
                }
                if (!coarse_residual.empty()) {
                    Real cx = 0.0, cy = 0.0;
                    for (int id : block) {
                        cx += db.nodes[id].x;
                        cy += db.nodes[id].y;
                    }
                    cx /= static_cast<Real>(block.size());
                    cy /= static_cast<Real>(block.size());
                    std::vector<std::pair<Real, int>> destinations;
                    destinations.reserve(coarse_residual.size());
                    for (int bin = 0; bin < static_cast<int>(coarse_residual.size());
                         ++bin) {
                        if (coarse_residual[bin] <= 0.0) continue;
                        const int bx = bin % coarse_x;
                        const int by = bin / coarse_x;
                        const Real tx = db.xl + (bx + 0.5) *
                            (db.xh - db.xl) / coarse_x;
                        const Real ty = db.yl + (by + 0.5) *
                            (db.yh - db.yl) / coarse_y;
                        destinations.emplace_back(std::hypot(tx - cx, ty - cy), bin);
                    }
                    std::sort(destinations.begin(), destinations.end());
                    const std::size_t destination_count =
                        std::min<std::size_t>(8, destinations.size());
                    for (std::size_t di = 0; di < destination_count; ++di) {
                        const int bin = destinations[di].second;
                        const int bx = bin % coarse_x;
                        const int by = bin / coarse_x;
                        const Real tx = db.xl + (bx + 0.5) *
                            (db.xh - db.xl) / coarse_x;
                        const Real ty = db.yl + (by + 0.5) *
                            (db.yh - db.yl) / coarse_y;
                        const Real dx = std::clamp(
                            tx - cx,
                            -config.step_bins * density.bin_width(),
                            config.step_bins * density.bin_width());
                        const Real dy = std::clamp(
                            ty - cy,
                            -config.step_bins * density.bin_height(),
                            config.step_bins * density.bin_height());
                        for (int line = 0; line < config.line_search_steps; ++line) {
                            const Real scale = std::ldexp(1.0, -line);
                            consider_block(dx * scale, dy * scale);
                            if (config.axis_separated) {
                                consider_block(dx * scale, 0.0);
                                consider_block(0.0, dy * scale);
                            }
                        }
                    }
                }
                if (config.breakpoint_oracle) {
                    for (int id : block) {
                        const Node& node = db.nodes[id];
                        for (Real x : nearest_breakpoints(
                                 node.x, node.width, db.xl, db.xh,
                                 density.bin_width(), 1.0))
                            consider_block(x - node.x, 0.0);
                        for (Real x : nearest_breakpoints(
                                 node.x, node.width, db.xl, db.xh,
                                 density.bin_width(), -1.0))
                            consider_block(x - node.x, 0.0);
                        for (Real y : nearest_breakpoints(
                                 node.y, node.height, db.yl, db.yh,
                                 density.bin_height(), 1.0))
                            consider_block(0.0, y - node.y);
                        for (Real y : nearest_breakpoints(
                                 node.y, node.height, db.yl, db.yh,
                                 density.bin_height(), -1.0))
                            consider_block(0.0, y - node.y);
                    }
                }
                if (!found) continue;
                std::vector<DensityNodeMove> moves;
                moves.reserve(block.size());
                for (int id : block) {
                    const Node& node = db.nodes[id];
                    moves.push_back({id,
                        std::clamp(node.x + best_dx,
                                   db.xl + 0.5 * node.width,
                                   db.xh - 0.5 * node.width),
                        std::clamp(node.y + best_dy,
                                   db.yl + 0.5 * node.height,
                                   db.yh - 0.5 * node.height)});
                }
                DensityMove committed = density.evaluate_group_move(moves);
                density.commit_move(committed);
                if (!coarse_residual.empty()) {
                    for (const auto& [fine, delta] : committed.area_changes) {
                        const int fx = fine % density.bins_x();
                        const int fy = fine / density.bins_x();
                        const int cx = fx * coarse_x / density.bins_x();
                        const int cy = fy * coarse_y / density.bins_y();
                        coarse_residual[cy * coarse_x + cx] -= delta;
                    }
                }
                for (const DensityNodeMove& move : moves) {
                    db.nodes[move.node_id].x = move.x;
                    db.nodes[move.node_id].y = move.y;
                    block_touched[move.node_id] = 1;
                }
                current_hpwl += best_hpwl_delta;
                ++sweep_moves;
                ++stats.moves;
            }
        }
        if (config.node_moves) {
        for (int id : db.movable_ids) {
            Node& node = db.nodes[id];
            bool found = false;
            Real best_x = node.x, best_y = node.y;
            Real best_overflow_delta = 0.0;
            Real best_hpwl_delta = std::numeric_limits<Real>::infinity();

            auto consider = [&](Real nx, Real ny) {
                if (std::abs(nx - node.x) <= 1.0e-10 &&
                    std::abs(ny - node.y) <= 1.0e-10) return;
                ++stats.candidates;
                DensityMove move = density.evaluate_move(id, nx, ny);
                if (move.overflow_area_delta >= -tolerance) return;
                const Real hd = affected_hpwl_delta(db, id, nx, ny, incident[id]);
                if (current_hpwl + hd > budget + 1.0e-6) return;
                if (!found || move.overflow_area_delta < best_overflow_delta - tolerance ||
                    (std::abs(move.overflow_area_delta - best_overflow_delta) <= tolerance &&
                     hd < best_hpwl_delta)) {
                    found = true;
                    best_x = nx;
                    best_y = ny;
                    best_overflow_delta = move.overflow_area_delta;
                    best_hpwl_delta = hd;
                }
            };

            for (int line = 0; line < config.line_search_steps; ++line) {
                const Real scale = std::ldexp(1.0, -line);
                const Real dx = config.step_bins * density.bin_width() * scale;
                const Real dy = config.step_bins * density.bin_height() * scale;
                const Real xs[2] = {
                    std::clamp(node.x - dx, db.xl + 0.5 * node.width,
                               db.xh - 0.5 * node.width),
                    std::clamp(node.x + dx, db.xl + 0.5 * node.width,
                               db.xh - 0.5 * node.width)};
                const Real ys[2] = {
                    std::clamp(node.y - dy, db.yl + 0.5 * node.height,
                               db.yh - 0.5 * node.height),
                    std::clamp(node.y + dy, db.yl + 0.5 * node.height,
                               db.yh - 0.5 * node.height)};
                for (Real x : xs) consider(x, node.y);
                for (Real y : ys) consider(node.x, y);
                if (!config.axis_separated) {
                    for (Real x : xs) for (Real y : ys) consider(x, y);
                }
            }
            if (config.breakpoint_oracle) {
                const auto xp = nearest_breakpoints(
                    node.x, node.width, db.xl, db.xh, density.bin_width(), 1.0);
                const auto xm = nearest_breakpoints(
                    node.x, node.width, db.xl, db.xh, density.bin_width(), -1.0);
                const auto yp = nearest_breakpoints(
                    node.y, node.height, db.yl, db.yh, density.bin_height(), 1.0);
                const auto ym = nearest_breakpoints(
                    node.y, node.height, db.yl, db.yh, density.bin_height(), -1.0);
                for (Real x : xp) consider(x, node.y);
                for (Real x : xm) consider(x, node.y);
                for (Real y : yp) consider(node.x, y);
                for (Real y : ym) consider(node.x, y);
            }
            if (!found) continue;
            DensityMove committed = density.evaluate_move(id, best_x, best_y);
            density.commit_move(committed);
            node.x = best_x;
            node.y = best_y;
            current_hpwl += best_hpwl_delta;
            ++sweep_moves;
            ++stats.moves;
        }
        }
        current = density.evaluate(0.0, 1.0, nullptr, nullptr);
        current_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        stats.trajectory.push_back({sweep, sweep_moves, current_hpwl, current.overflow});
        ++stats.sweeps;
        if (sweep_moves == 0) break;
    }
    stats.final_hpwl = current_hpwl;
    stats.final_overflow = current.overflow;
    stats.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return stats;
}

}  // namespace ea
