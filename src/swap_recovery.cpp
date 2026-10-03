#include "epsilon_active/swap_recovery.hpp"

#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <queue>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ea {
namespace {

struct BucketKey {
    std::int64_t width = 0;
    std::int64_t height = 0;
    int bin = 0;

    bool operator==(const BucketKey& other) const noexcept {
        return width == other.width && height == other.height && bin == other.bin;
    }
};

struct BucketHash {
    std::size_t operator()(const BucketKey& key) const noexcept {
        std::size_t value = std::hash<std::int64_t>{}(key.width);
        value ^= std::hash<std::int64_t>{}(key.height) + 0x9e3779b9U +
                 (value << 6U) + (value >> 2U);
        value ^= std::hash<int>{}(key.bin) + 0x9e3779b9U +
                 (value << 6U) + (value >> 2U);
        return value;
    }
};

std::int64_t dimension_key(Real value) {
    return static_cast<std::int64_t>(std::llround(value * 1000000.0));
}

int spatial_bin(const Database& db, const ExactOverlapDensity& density,
                Real x, Real y) {
    const int bx = std::clamp(
        static_cast<int>(std::floor((x - db.xl) / density.bin_width())),
        0, density.bins_x() - 1);
    const int by = std::clamp(
        static_cast<int>(std::floor((y - db.yl) / density.bin_height())),
        0, density.bins_y() - 1);
    return by * density.bins_x() + bx;
}

struct AreaChange { int bin = -1; Real delta = 0.0; };

void rectangle_changes(const Database& db, const ExactOverlapDensity& density,
                       const Node& node, Real x, Real y, Real sign,
                       std::vector<AreaChange>& out) {
    const Real left = x - 0.5 * node.width;
    const Real right = x + 0.5 * node.width;
    const Real bottom = y - 0.5 * node.height;
    const Real top = y + 0.5 * node.height;
    const int x0 = std::clamp(static_cast<int>(std::floor((left - db.xl) /
        density.bin_width())), 0, density.bins_x() - 1);
    const int x1 = std::clamp(static_cast<int>(std::floor((std::nextafter(
        right, left) - db.xl) / density.bin_width())), 0, density.bins_x() - 1);
    const int y0 = std::clamp(static_cast<int>(std::floor((bottom - db.yl) /
        density.bin_height())), 0, density.bins_y() - 1);
    const int y1 = std::clamp(static_cast<int>(std::floor((std::nextafter(
        top, bottom) - db.yl) / density.bin_height())), 0, density.bins_y() - 1);
    for (int by = y0; by <= y1; ++by) {
        const Real bin_bottom = db.yl + by * density.bin_height();
        const Real oy = std::max<Real>(0.0, std::min(top, bin_bottom + density.bin_height()) -
            std::max(bottom, bin_bottom));
        for (int bx = x0; bx <= x1; ++bx) {
            const Real bin_left = db.xl + bx * density.bin_width();
            const Real ox = std::max<Real>(0.0, std::min(right, bin_left + density.bin_width()) -
                std::max(left, bin_left));
            if (ox > 0.0 && oy > 0.0) out.push_back({by * density.bins_x() + bx,
                                                       sign * ox * oy});
        }
    }
}

Real swap_overflow_delta(const Database& db, const ExactOverlapDensity& density,
                         const Node& a, const Node& b,
                         const std::vector<Real>& occupancy,
                         std::vector<AreaChange>& changes) {
    changes.clear();
    rectangle_changes(db, density, a, a.x, a.y, -1.0, changes);
    rectangle_changes(db, density, b, b.x, b.y, -1.0, changes);
    rectangle_changes(db, density, a, b.x, b.y, 1.0, changes);
    rectangle_changes(db, density, b, a.x, a.y, 1.0, changes);
    std::sort(changes.begin(), changes.end(), [](const AreaChange& x, const AreaChange& y) {
        return x.bin < y.bin;
    });
    Real result = 0.0;
    const Real capacity = density.target_density() * density.bin_area();
    std::size_t write = 0;
    for (std::size_t i = 0; i < changes.size();) {
        const int bin = changes[i].bin;
        Real delta = 0.0;
        while (i < changes.size() && changes[i].bin == bin) delta += changes[i++].delta;
        const Real old_excess = std::max<Real>(occupancy[bin] - capacity, 0.0);
        const Real new_excess = std::max<Real>(occupancy[bin] + delta - capacity, 0.0);
        result += new_excess - old_excess;
        if (std::abs(delta) > 1.0e-14) changes[write++] = {bin, delta};
    }
    changes.resize(write);
    return result;
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

std::vector<std::pair<Real, Real>> net_targets(const Database& db) {
    std::vector<Real> sum_x(db.nodes.size(), 0.0);
    std::vector<Real> sum_y(db.nodes.size(), 0.0);
    std::vector<int> count(db.nodes.size(), 0);
    for (const Net& net : db.nets) {
        if (net.pin_count == 0) continue;
        Real cx = 0.0;
        Real cy = 0.0;
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            const Pin& pin = db.pins[p];
            cx += db.nodes[pin.node].x + pin.offset_x;
            cy += db.nodes[pin.node].y + pin.offset_y;
        }
        cx /= static_cast<Real>(net.pin_count);
        cy /= static_cast<Real>(net.pin_count);
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            const int id = db.pins[p].node;
            sum_x[id] += cx;
            sum_y[id] += cy;
            ++count[id];
        }
    }
    std::vector<std::pair<Real, Real>> result(db.nodes.size());
    for (std::size_t id = 0; id < db.nodes.size(); ++id) {
        result[id] = count[id] > 0
            ? std::pair<Real, Real>{sum_x[id] / count[id], sum_y[id] / count[id]}
            : std::pair<Real, Real>{db.nodes[id].x, db.nodes[id].y};
    }
    return result;
}

Real swap_hpwl_delta(const Database& db, int a, int b,
                     const std::vector<std::vector<int>>& node_nets,
                     int degree_limit) {
    std::vector<int> nets = node_nets[a];
    nets.insert(nets.end(), node_nets[b].begin(), node_nets[b].end());
    std::sort(nets.begin(), nets.end());
    nets.erase(std::unique(nets.begin(), nets.end()), nets.end());
    Real delta = 0.0;
    for (int net_id : nets) {
        const Net& net = db.nets[net_id];
        if (net.pin_count < 2 || net.pin_count > static_cast<std::size_t>(degree_limit)) {
            continue;
        }
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
            Real candidate_x = node.x;
            Real candidate_y = node.y;
            if (pin.node == a) {
                candidate_x = db.nodes[b].x;
                candidate_y = db.nodes[b].y;
            } else if (pin.node == b) {
                candidate_x = db.nodes[a].x;
                candidate_y = db.nodes[a].y;
            }
            old_min_x = std::min(old_min_x, node.x + pin.offset_x);
            old_max_x = std::max(old_max_x, node.x + pin.offset_x);
            old_min_y = std::min(old_min_y, node.y + pin.offset_y);
            old_max_y = std::max(old_max_y, node.y + pin.offset_y);
            new_min_x = std::min(new_min_x, candidate_x + pin.offset_x);
            new_max_x = std::max(new_max_x, candidate_x + pin.offset_x);
            new_min_y = std::min(new_min_y, candidate_y + pin.offset_y);
            new_max_y = std::max(new_max_y, candidate_y + pin.offset_y);
        }
        delta += net.weight * ((new_max_x - new_min_x) + (new_max_y - new_min_y) -
                               (old_max_x - old_min_x) - (old_max_y - old_min_y));
    }
    return delta;
}

Real node_hpwl_at(const Database& db, int node_id, Real new_x, Real new_y,
                  const std::vector<std::vector<int>>& node_nets) {
    Real total = 0.0;
    for (int net_id : node_nets[node_id]) {
        const Net& net = db.nets[net_id];
        if (net.pin_count < 2) continue;
        Real min_x = std::numeric_limits<Real>::infinity();
        Real max_x = -std::numeric_limits<Real>::infinity();
        Real min_y = std::numeric_limits<Real>::infinity();
        Real max_y = -std::numeric_limits<Real>::infinity();
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            const Pin& pin = db.pins[p];
            const Node& node = db.nodes[pin.node];
            const Real x = (pin.node == node_id ? new_x : node.x) + pin.offset_x;
            const Real y = (pin.node == node_id ? new_y : node.y) + pin.offset_y;
            min_x = std::min(min_x, x);
            max_x = std::max(max_x, x);
            min_y = std::min(min_y, y);
            max_y = std::max(max_y, y);
        }
        total += net.weight * ((max_x - min_x) + (max_y - min_y));
    }
    return total;
}

std::vector<std::vector<int>> net_disjoint_batches(
    const std::vector<int>& nodes,
    const std::vector<std::vector<int>>& node_nets) {
    std::vector<int> order = nodes;
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (node_nets[a].size() != node_nets[b].size()) {
            return node_nets[a].size() > node_nets[b].size();
        }
        return a < b;
    });
    std::unordered_map<int, std::vector<int>> net_colors;
    std::vector<std::vector<int>> batches;
    std::vector<int> mark;
    int stamp = 0;
    for (int id : order) {
        ++stamp;
        if (mark.size() < batches.size() + 1) mark.resize(batches.size() + 1, 0);
        for (int net_id : node_nets[id]) {
            const auto found = net_colors.find(net_id);
            if (found == net_colors.end()) continue;
            for (int color : found->second) mark[color] = stamp;
        }
        int color = 0;
        while (color < static_cast<int>(batches.size()) && mark[color] == stamp) {
            ++color;
        }
        if (color == static_cast<int>(batches.size())) {
            batches.emplace_back();
            mark.push_back(0);
        }
        batches[color].push_back(id);
        for (int net_id : node_nets[id]) net_colors[net_id].push_back(color);
    }
    return batches;
}

struct AssignmentResult {
    bool complete = false;
    int moved = 0;
    std::int64_t bids = 0;
    Real predicted_delta = 0.0;
    std::vector<int> anchor_of;
};

AssignmentResult assign_exact_anchors(
    const Database& db, const ExactOverlapDensity& density,
    const std::vector<int>& batch,
    const std::vector<std::vector<int>>& node_nets,
    const std::vector<std::pair<Real, Real>>& targets,
    const SwapRecoveryConfig& config) {
    AssignmentResult result;
    const int count = static_cast<int>(batch.size());
    if (count < 2) return result;
    std::unordered_map<int, std::vector<int>> bins;
    bins.reserve(batch.size() * 2);
    for (int anchor = 0; anchor < count; ++anchor) {
        const Node& node = db.nodes[batch[anchor]];
        bins[spatial_bin(db, density, node.x, node.y)].push_back(anchor);
    }

    std::vector<std::vector<std::pair<int, Real>>> edges(count);
    const auto& readonly_bins = bins;
    #pragma omp parallel for schedule(dynamic, 64) if(count > 256)
    for (int i = 0; i < count; ++i) {
        const int id = batch[i];
        const Node& node = db.nodes[id];
        const Real own_cost = node_hpwl_at(db, id, node.x, node.y, node_nets);
        std::vector<int> candidates = {i};
        const int target_bin = spatial_bin(
            db, density, targets[id].first, targets[id].second);
        const int target_x = target_bin % density.bins_x();
        const int target_y = target_bin / density.bins_x();
        for (int radius = 0;
             radius <= config.assignment_radius_bins &&
             static_cast<int>(candidates.size()) < config.assignment_candidates;
             ++radius) {
            const int x0 = std::max(0, target_x - radius);
            const int x1 = std::min(density.bins_x() - 1, target_x + radius);
            const int y0 = std::max(0, target_y - radius);
            const int y1 = std::min(density.bins_y() - 1, target_y + radius);
            for (int by = y0; by <= y1 &&
                 static_cast<int>(candidates.size()) < config.assignment_candidates; ++by) {
                for (int bx = x0; bx <= x1 &&
                     static_cast<int>(candidates.size()) < config.assignment_candidates; ++bx) {
                    if (radius > 0 && bx > x0 && bx < x1 && by > y0 && by < y1) {
                        continue;
                    }
                    const auto found = readonly_bins.find(by * density.bins_x() + bx);
                    if (found == readonly_bins.end()) continue;
                    for (int anchor : found->second) {
                        if (std::find(candidates.begin(), candidates.end(), anchor) ==
                            candidates.end()) {
                            candidates.push_back(anchor);
                        }
                        if (static_cast<int>(candidates.size()) >=
                            config.assignment_candidates) break;
                    }
                }
            }
        }
        edges[i].reserve(candidates.size());
        for (int anchor : candidates) {
            const Node& position = db.nodes[batch[anchor]];
            const Real cost = node_hpwl_at(
                db, id, position.x, position.y, node_nets) - own_cost;
            edges[i].emplace_back(anchor, cost);
        }
    }

    std::vector<Real> prices(count, 0.0);
    std::vector<int> owner(count, -1);
    result.anchor_of.assign(count, -1);
    std::deque<int> queue;
    for (int i = 0; i < count; ++i) queue.push_back(i);
    const std::int64_t maximum_bids = static_cast<std::int64_t>(
        config.assignment_max_bids_per_node) * count;
    while (!queue.empty() && result.bids < maximum_bids) {
        const int person = queue.front();
        queue.pop_front();
        int best_anchor = -1;
        Real best_value = -std::numeric_limits<Real>::infinity();
        Real second_value = -std::numeric_limits<Real>::infinity();
        for (const auto& [anchor, cost] : edges[person]) {
            const Real value = -cost - prices[anchor];
            if (value > best_value) {
                second_value = best_value;
                best_value = value;
                best_anchor = anchor;
            } else if (value > second_value) {
                second_value = value;
            }
        }
        if (best_anchor < 0) break;
        if (!std::isfinite(second_value)) {
            second_value = best_value - config.assignment_epsilon;
        }
        prices[best_anchor] += best_value - second_value +
                               config.assignment_epsilon;
        const int displaced = owner[best_anchor];
        owner[best_anchor] = person;
        result.anchor_of[person] = best_anchor;
        if (displaced >= 0) {
            result.anchor_of[displaced] = -1;
            queue.push_back(displaced);
        }
        ++result.bids;
    }
    if (!queue.empty()) return result;
    result.complete = true;
    for (int i = 0; i < count; ++i) {
        if (result.anchor_of[i] < 0) {
            result.complete = false;
            return result;
        }
        if (result.anchor_of[i] != i) ++result.moved;
        const auto found = std::find_if(
            edges[i].begin(), edges[i].end(), [&](const auto& edge) {
                return edge.first == result.anchor_of[i];
            });
        if (found == edges[i].end()) {
            result.complete = false;
            return result;
        }
        result.predicted_delta += found->second;
    }
    return result;
}

std::uint64_t hilbert_index(int n, int x, int y) {
    std::uint64_t distance = 0;
    for (int scale = n / 2; scale > 0; scale /= 2) {
        const int rx = (x & scale) != 0;
        const int ry = (y & scale) != 0;
        distance += static_cast<std::uint64_t>(scale) * scale * ((3 * rx) ^ ry);
        if (ry == 0) {
            if (rx == 1) {
                x = n - 1 - x;
                y = n - 1 - y;
            }
            std::swap(x, y);
        }
    }
    return distance;
}

std::uint64_t position_hilbert(const Database& db,
                               const ExactOverlapDensity& density,
                               Real x, Real y) {
    const int bin = spatial_bin(db, density, x, y);
    return hilbert_index(density.bins_x(), bin % density.bins_x(),
                         bin / density.bins_x());
}

}  // namespace

SwapRecoveryStats recover_hpwl_with_equal_shape_swaps(
    Database& db, ExactOverlapDensity& density,
    const SwapRecoveryConfig& config) {
    if (config.sweeps < 0 || config.permutation_sweeps < 0 ||
        config.assignment_sweeps < 0 || config.assignment_radius_bins < 0 ||
        config.assignment_candidates < 2 ||
        config.assignment_max_bids_per_node <= 0 ||
        config.assignment_epsilon <= 0.0 ||
        config.radius_bins < 0 || config.candidates <= 0 ||
        config.exact_shortlist < 0 || config.degree_limit < 2) {
        throw std::invalid_argument("invalid equal-shape swap recovery configuration");
    }
    SwapRecoveryStats stats;
    if (config.sweeps == 0 && config.permutation_sweeps == 0 &&
        config.assignment_sweeps == 0) return stats;
    if (config.permutation_sweeps > 0 &&
        (density.bins_x() != density.bins_y() ||
         (density.bins_x() & (density.bins_x() - 1)) != 0)) {
        throw std::invalid_argument(
            "shape permutation requires a square power-of-two density grid");
    }
    const auto started = std::chrono::steady_clock::now();
    ExactHpwl hpwl(db);
    DensityMetrics density_metrics = density.evaluate(0.0, 1.0, nullptr, nullptr);
    stats.initial_overflow = density_metrics.overflow;
    stats.initial_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    stats.objective_evaluations += 2;
    const auto node_nets = build_node_nets(db);
    std::vector<Real> occupancy = density.occupancy();
    const Real capacity_tolerance = 1.0e-10 * density.bin_area();
    std::vector<unsigned char> swapped(db.nodes.size(), 0);

    for (int sweep = 0; sweep < config.sweeps; ++sweep) {
        std::fill(swapped.begin(), swapped.end(), 0);
        const auto targets = net_targets(db);
        std::unordered_map<BucketKey, std::vector<int>, BucketHash> buckets;
        buckets.reserve(db.movable_ids.size() * 2);
        std::unordered_map<int, std::vector<int>> spatial_buckets;
        if (config.density_guided && config.allow_unequal_shapes) {
            spatial_buckets.reserve(db.movable_ids.size() / 4);
        }
        for (int id : db.movable_ids) {
            const Node& node = db.nodes[id];
            buckets[{dimension_key(node.width), dimension_key(node.height),
                     spatial_bin(db, density, node.x, node.y)}].push_back(id);
            if (config.density_guided && config.allow_unequal_shapes) {
                spatial_buckets[spatial_bin(db, density, node.x, node.y)].push_back(id);
            }
        }
        std::vector<int> order = db.movable_ids;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return db.node_pin_count[a] > db.node_pin_count[b];
        });
        int sweep_swaps = 0;
        const Real sweep_start_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        std::vector<int> underfull_bins;
        std::vector<int> nearest_underfull;
        if (config.density_guided) {
            const Real cap = density.target_density() * density.bin_area();
            const auto& occ = density.occupancy();
            for (int bin = 0; bin < static_cast<int>(occ.size()); ++bin) {
                if (occ[bin] < cap - capacity_tolerance) underfull_bins.push_back(bin);
            }
            nearest_underfull.assign(occ.size(), -1);
            std::queue<int> queue;
            for (int bin : underfull_bins) {
                nearest_underfull[bin] = bin;
                queue.push(bin);
            }
            while (!queue.empty()) {
                const int bin = queue.front();
                queue.pop();
                const int bx = bin % density.bins_x();
                const int by = bin / density.bins_x();
                const int neighbors[4] = {
                    bx > 0 ? bin - 1 : -1,
                    bx + 1 < density.bins_x() ? bin + 1 : -1,
                    by > 0 ? bin - density.bins_x() : -1,
                    by + 1 < density.bins_y() ? bin + density.bins_x() : -1};
                for (int next : neighbors) {
                    if (next >= 0 && nearest_underfull[next] < 0) {
                        nearest_underfull[next] = nearest_underfull[bin];
                        queue.push(next);
                    }
                }
            }
        }
        for (int a : order) {
            if (swapped[a]) continue;
            const Node& node_a = db.nodes[a];
            int target_bin = spatial_bin(
                db, density, targets[a].first, targets[a].second);
            if (config.density_guided && !underfull_bins.empty()) {
                const int current_bin = spatial_bin(db, density, node_a.x, node_a.y);
                if (nearest_underfull[current_bin] >= 0) target_bin = nearest_underfull[current_bin];
            }
            const int target_x = target_bin % density.bins_x();
            const int target_y = target_bin / density.bins_x();
            bool found = false;
            int best_b = -1;
            Real best_delta = 0.0;
            Real best_overflow_delta = 0.0;
            Real best_score = std::numeric_limits<Real>::infinity();
            std::vector<AreaChange> best_changes;
            int tested = 0;
            std::vector<std::pair<Real, int>> pool;
            pool.reserve(config.candidates);
            std::unordered_set<int> seen_candidates;
            seen_candidates.reserve(config.candidates * 2 + 8);
            if (config.net_aware) {
                for (int net_id : node_nets[a]) {
                    const Net& net = db.nets[net_id];
                    if (net.pin_count > static_cast<std::size_t>(config.degree_limit)) {
                        continue;
                    }
                    for (std::size_t p = net.pin_begin;
                         p < net.pin_begin + net.pin_count &&
                         tested < config.candidates; ++p) {
                        const int b = db.pins[p].node;
                        if (b == a || swapped[b] ||
                            !seen_candidates.insert(b).second) continue;
                        const Node& node_b = db.nodes[b];
                        if (node_a.width != node_b.width ||
                            node_a.height != node_b.height) continue;
                        ++tested;
                        pool.emplace_back(0.0, b);
                    }
                    if (tested >= config.candidates) break;
                }
            }
            for (int radius = 0; radius <= config.radius_bins && tested < config.candidates; ++radius) {
                const int x0 = std::max(0, target_x - radius);
                const int x1 = std::min(density.bins_x() - 1, target_x + radius);
                const int y0 = std::max(0, target_y - radius);
                const int y1 = std::min(density.bins_y() - 1, target_y + radius);
                for (int by = y0; by <= y1 && tested < config.candidates; ++by) {
                    for (int bx = x0; bx <= x1 && tested < config.candidates; ++bx) {
                        if (radius > 0 && bx > x0 && bx < x1 && by > y0 && by < y1) continue;
                        const int candidate_bin = by * density.bins_x() + bx;
                        std::vector<int> candidates;
                        if (config.density_guided && config.allow_unequal_shapes) {
                            const auto bucket = spatial_buckets.find(candidate_bin);
                            if (bucket != spatial_buckets.end()) candidates = bucket->second;
                        } else {
                            const BucketKey key{dimension_key(node_a.width),
                                                dimension_key(node_a.height), candidate_bin};
                            const auto bucket = buckets.find(key);
                            if (bucket != buckets.end()) candidates = bucket->second;
                        }
                        for (int b : candidates) {
                            if (tested >= config.candidates) break;
                            if (a == b || swapped[b] ||
                                !seen_candidates.insert(b).second) continue;
                            const Node& node_b = db.nodes[b];
                            if ((!config.density_guided || !config.allow_unequal_shapes) &&
                                (node_a.width != node_b.width || node_a.height != node_b.height)) {
                                continue;
                            }
                            if (spatial_bin(db, density, node_b.x, node_b.y) != candidate_bin) {
                                continue;
                            }
                            if (b == a ||
                                node_b.x < db.xl + 0.5 * node_a.width ||
                                node_b.x > db.xh - 0.5 * node_a.width ||
                                node_b.y < db.yl + 0.5 * node_a.height ||
                                node_b.y > db.yh - 0.5 * node_a.height ||
                                node_a.x < db.xl + 0.5 * node_b.width ||
                                node_a.x > db.xh - 0.5 * node_b.width ||
                                node_a.y < db.yl + 0.5 * node_b.height ||
                                node_a.y > db.yh - 0.5 * node_b.height) {
                                continue;
                            }
                            ++tested;
                            const Real a_to_b_x = node_b.x - targets[a].first;
                            const Real a_to_b_y = node_b.y - targets[a].second;
                            const Real b_to_a_x = node_a.x - targets[b].first;
                            const Real b_to_a_y = node_a.y - targets[b].second;
                            pool.emplace_back(
                                a_to_b_x * a_to_b_x + a_to_b_y * a_to_b_y +
                                b_to_a_x * b_to_a_x + b_to_a_y * b_to_a_y,
                                b);
                        }
                    }
                }
            }
            std::sort(pool.begin(), pool.end());
            if (config.exact_shortlist > 0 &&
                static_cast<int>(pool.size()) > config.exact_shortlist) {
                pool.resize(config.exact_shortlist);
            }
            for (const auto& [proxy, b] : pool) {
                (void)proxy;
                const Real delta = swap_hpwl_delta(
                    db, a, b, node_nets, config.degree_limit);
                if (config.density_guided) {
                    std::vector<AreaChange> changes;
                    const Real overflow_delta = swap_overflow_delta(
                        db, density, node_a, db.nodes[b], occupancy, changes);
                    if (overflow_delta >= -capacity_tolerance ||
                        delta > config.density_hpwl_budget * sweep_start_hpwl) continue;
                    const Real score = delta / std::max<Real>(sweep_start_hpwl, 1.0) +
                        config.density_score_weight * overflow_delta /
                        std::max<Real>(density.bin_area(), 1.0);
                    if (!found || score < best_score) {
                        found = true;
                        best_b = b;
                        best_delta = delta;
                        best_overflow_delta = overflow_delta;
                        best_score = score;
                        best_changes = std::move(changes);
                    }
                } else if (delta < best_delta - 1.0e-12) {
                    found = true;
                    best_b = b;
                    best_delta = delta;
                }
            }
            if (!found) continue;
            std::swap(db.nodes[a].x, db.nodes[best_b].x);
            std::swap(db.nodes[a].y, db.nodes[best_b].y);
            swapped[a] = 1;
            swapped[best_b] = 1;
            if (config.density_guided) {
                for (const AreaChange& change : best_changes) occupancy[change.bin] += change.delta;
                (void)best_overflow_delta;
            }
            ++sweep_swaps;
            ++stats.swaps;
        }
        density_metrics = density.evaluate(0.0, 1.0, nullptr, nullptr);
        const Real exact_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        stats.objective_evaluations += 2;
        stats.trajectory.push_back({sweep, sweep_swaps, exact_hpwl,
                                    density_metrics.overflow});
        ++stats.sweeps;
        if (sweep_swaps == 0) break;
    }
    Real current_hpwl = stats.trajectory.empty()
        ? stats.initial_hpwl : stats.trajectory.back().hpwl;
    for (int sweep = 0; sweep < config.permutation_sweeps; ++sweep) {
        const auto targets = net_targets(db);
        std::unordered_map<BucketKey, std::vector<int>, BucketHash> shapes;
        shapes.reserve(64);
        for (int id : db.movable_ids) {
            const Node& node = db.nodes[id];
            shapes[{dimension_key(node.width), dimension_key(node.height), 0}]
                .push_back(id);
        }
        std::vector<Real> old_x(db.movable_ids.size());
        std::vector<Real> old_y(db.movable_ids.size());
        for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
            const Node& node = db.nodes[db.movable_ids[i]];
            old_x[i] = node.x;
            old_y[i] = node.y;
        }
        int moved_nodes = 0;
        for (auto& [shape, nodes] : shapes) {
            (void)shape;
            if (nodes.size() < 2) continue;
            std::vector<std::pair<Real, Real>> positions;
            positions.reserve(nodes.size());
            for (int id : nodes) positions.emplace_back(db.nodes[id].x, db.nodes[id].y);
            std::sort(nodes.begin(), nodes.end(), [&](int a, int b) {
                const std::uint64_t ha = position_hilbert(
                    db, density, targets[a].first, targets[a].second);
                const std::uint64_t hb = position_hilbert(
                    db, density, targets[b].first, targets[b].second);
                return ha == hb ? a < b : ha < hb;
            });
            std::sort(positions.begin(), positions.end(), [&](const auto& a, const auto& b) {
                const std::uint64_t ha = position_hilbert(db, density, a.first, a.second);
                const std::uint64_t hb = position_hilbert(db, density, b.first, b.second);
                if (ha != hb) return ha < hb;
                if (a.first != b.first) return a.first < b.first;
                return a.second < b.second;
            });
            for (std::size_t i = 0; i < nodes.size(); ++i) {
                Node& node = db.nodes[nodes[i]];
                if (node.x != positions[i].first || node.y != positions[i].second) {
                    ++moved_nodes;
                }
                node.x = positions[i].first;
                node.y = positions[i].second;
            }
        }
        const Real candidate_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        ++stats.objective_evaluations;
        if (candidate_hpwl >= current_hpwl - 1.0e-12) {
            for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
                Node& node = db.nodes[db.movable_ids[i]];
                node.x = old_x[i];
                node.y = old_y[i];
            }
            break;
        }
        density_metrics = density.evaluate(0.0, 1.0, nullptr, nullptr);
        ++stats.objective_evaluations;
        current_hpwl = candidate_hpwl;
        ++stats.permutation_sweeps;
        stats.permuted_nodes += moved_nodes;
        stats.trajectory.push_back({stats.sweeps, moved_nodes, candidate_hpwl,
                                    density_metrics.overflow});
        ++stats.sweeps;
    }
    for (int sweep = 0; sweep < config.assignment_sweeps; ++sweep) {
        const auto targets = net_targets(db);
        std::unordered_map<BucketKey, std::vector<int>, BucketHash> shape_map;
        shape_map.reserve(64);
        for (int id : db.movable_ids) {
            const Node& node = db.nodes[id];
            shape_map[{dimension_key(node.width), dimension_key(node.height), 0}]
                .push_back(id);
        }
        std::vector<std::pair<BucketKey, std::vector<int>*>> shapes;
        shapes.reserve(shape_map.size());
        for (auto& [shape, nodes] : shape_map) shapes.emplace_back(shape, &nodes);
        std::sort(shapes.begin(), shapes.end(), [](const auto& a, const auto& b) {
            if (a.first.width != b.first.width) {
                return a.first.width < b.first.width;
            }
            return a.first.height < b.first.height;
        });

        int sweep_moved = 0;
        for (auto& [shape, shape_nodes] : shapes) {
            (void)shape;
            if (shape_nodes->size() < 2) continue;
            const auto batches = net_disjoint_batches(*shape_nodes, node_nets);
            for (const std::vector<int>& batch : batches) {
                if (batch.size() < 2) continue;
                AssignmentResult assignment = assign_exact_anchors(
                    db, density, batch, node_nets, targets, config);
                stats.assignment_bids += assignment.bids;
                if (!assignment.complete || assignment.moved == 0 ||
                    assignment.predicted_delta >= -1.0e-12) {
                    continue;
                }
                std::vector<Real> old_x(batch.size());
                std::vector<Real> old_y(batch.size());
                for (std::size_t i = 0; i < batch.size(); ++i) {
                    old_x[i] = db.nodes[batch[i]].x;
                    old_y[i] = db.nodes[batch[i]].y;
                }
                for (std::size_t i = 0; i < batch.size(); ++i) {
                    Node& node = db.nodes[batch[i]];
                    const int anchor = assignment.anchor_of[i];
                    node.x = old_x[anchor];
                    node.y = old_y[anchor];
                }
                const Real candidate_hpwl = hpwl.evaluate(
                    0.0, 1.0, -1, nullptr, nullptr);
                const DensityMetrics candidate_density = density.evaluate(
                    0.0, 1.0, nullptr, nullptr);
                stats.objective_evaluations += 2;
                if (candidate_hpwl >= current_hpwl - 1.0e-12 ||
                    std::abs(candidate_density.overflow - density_metrics.overflow) >
                        1.0e-10) {
                    for (std::size_t i = 0; i < batch.size(); ++i) {
                        Node& node = db.nodes[batch[i]];
                        node.x = old_x[i];
                        node.y = old_y[i];
                    }
                    continue;
                }
                current_hpwl = candidate_hpwl;
                density_metrics = candidate_density;
                sweep_moved += assignment.moved;
                stats.assigned_nodes += assignment.moved;
                ++stats.assignment_batches;
            }
        }
        ++stats.assignment_sweeps;
        ++stats.sweeps;
        stats.trajectory.push_back({stats.sweeps - 1, sweep_moved, current_hpwl,
                                    density_metrics.overflow});
        if (sweep_moved == 0) break;
    }
    stats.final_hpwl = stats.trajectory.empty()
        ? stats.initial_hpwl : stats.trajectory.back().hpwl;
    stats.final_overflow = stats.trajectory.empty()
        ? stats.initial_overflow : stats.trajectory.back().overflow;
    stats.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return stats;
}

}  // namespace ea
