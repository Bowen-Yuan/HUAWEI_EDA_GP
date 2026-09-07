#include "epsilon_active/transport.hpp"

#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ea {
namespace {

struct AreaChange {
    int bin = -1;
    Real delta = 0.0;
};

struct MoveDelta {
    Real overflow_area = 0.0;
    Real energy = 0.0;
    std::vector<AreaChange> changes;
};

int lower_bin(Real coordinate, Real origin, Real pitch, int count) {
    return std::clamp(static_cast<int>(std::floor((coordinate - origin) / pitch)),
                      0, count - 1);
}

void rectangle_changes(const Database& db, const ExactOverlapDensity& density,
                       const Node& node, Real x, Real y, Real sign,
                       std::vector<AreaChange>& changes) {
    const Real left = x - 0.5 * node.width;
    const Real right = x + 0.5 * node.width;
    const Real bottom = y - 0.5 * node.height;
    const Real top = y + 0.5 * node.height;
    if (right <= db.xl || left >= db.xh || top <= db.yl || bottom >= db.yh) return;
    const int x0 = lower_bin(left, db.xl, density.bin_width(), density.bins_x());
    const int x1 = lower_bin(std::nextafter(right, left), db.xl,
                             density.bin_width(), density.bins_x());
    const int y0 = lower_bin(bottom, db.yl, density.bin_height(), density.bins_y());
    const int y1 = lower_bin(std::nextafter(top, bottom), db.yl,
                             density.bin_height(), density.bins_y());
    for (int by = y0; by <= y1; ++by) {
        const Real bin_bottom = db.yl + by * density.bin_height();
        const Real overlap_y = std::max<Real>(
            0.0, std::min(top, bin_bottom + density.bin_height()) -
                 std::max(bottom, bin_bottom));
        for (int bx = x0; bx <= x1; ++bx) {
            const Real bin_left = db.xl + bx * density.bin_width();
            const Real overlap_x = std::max<Real>(
                0.0, std::min(right, bin_left + density.bin_width()) -
                     std::max(left, bin_left));
            if (overlap_x > 0.0 && overlap_y > 0.0) {
                changes.push_back({by * density.bins_x() + bx,
                                   sign * overlap_x * overlap_y});
            }
        }
    }
}

MoveDelta evaluate_move(const Database& db, const ExactOverlapDensity& density,
                        const Node& node, Real new_x, Real new_y,
                        const std::vector<Real>& occupancy) {
    MoveDelta result;
    rectangle_changes(db, density, node, node.x, node.y, -1.0, result.changes);
    rectangle_changes(db, density, node, new_x, new_y, 1.0, result.changes);
    std::sort(result.changes.begin(), result.changes.end(),
              [](const AreaChange& a, const AreaChange& b) { return a.bin < b.bin; });
    std::size_t write = 0;
    for (std::size_t read = 0; read < result.changes.size();) {
        const int bin = result.changes[read].bin;
        Real delta = 0.0;
        while (read < result.changes.size() && result.changes[read].bin == bin) {
            delta += result.changes[read].delta;
            ++read;
        }
        if (std::abs(delta) > 1.0e-14) result.changes[write++] = {bin, delta};
    }
    result.changes.resize(write);

    const Real capacity = density.target_density() * density.bin_area();
    for (const AreaChange& change : result.changes) {
        const Real old_excess = std::max<Real>(occupancy[change.bin] - capacity, 0.0);
        const Real new_occupancy = occupancy[change.bin] + change.delta;
        const Real new_excess = std::max<Real>(new_occupancy - capacity, 0.0);
        result.overflow_area += new_excess - old_excess;
        result.energy += 0.5 / density.bin_area() *
            (new_excess * new_excess - old_excess * old_excess);
    }
    return result;
}

void apply_move(Node& node, Real x, Real y, const MoveDelta& delta,
                std::vector<Real>& occupancy) {
    for (const AreaChange& change : delta.changes) occupancy[change.bin] += change.delta;
    node.x = x;
    node.y = y;
}

enum class AtomicMoveResult { Accepted, Rejected, CapacityRejected };

AtomicMoveResult apply_atomic_group_move(
    Database& db, const ExactOverlapDensity& density,
    const std::vector<int>& members, Real dx, Real dy,
    Real maximum_overflow_delta, Real tolerance,
    bool require_componentwise_capacity, bool commit,
    std::vector<Real>& occupancy) {
    for (int id : members) {
        const Node& node = db.nodes[id];
        const Real x = node.x + dx;
        const Real y = node.y + dy;
        if (x < db.xl + 0.5 * node.width ||
            x > db.xh - 0.5 * node.width ||
            y < db.yl + 0.5 * node.height ||
            y > db.yh - 0.5 * node.height) {
            return AtomicMoveResult::Rejected;
        }
    }

    std::vector<std::pair<Real, Real>> old_positions;
    old_positions.reserve(members.size());
    std::vector<std::pair<int, Real>> saved_bins;
    Real overflow_delta = 0.0;
    for (int id : members) {
        Node& node = db.nodes[id];
        old_positions.emplace_back(node.x, node.y);
        const Real x = node.x + dx;
        const Real y = node.y + dy;
        MoveDelta delta = evaluate_move(db, density, node, x, y, occupancy);
        overflow_delta += delta.overflow_area;
        for (const AreaChange& change : delta.changes) {
            const auto saved = std::find_if(
                saved_bins.begin(), saved_bins.end(), [&](const auto& item) {
                    return item.first == change.bin;
                });
            if (saved == saved_bins.end()) {
                saved_bins.emplace_back(change.bin, occupancy[change.bin]);
            }
        }
        apply_move(node, x, y, delta, occupancy);
    }

    const bool aggregate_accepted = overflow_delta < -tolerance &&
        overflow_delta <= maximum_overflow_delta + tolerance;
    bool componentwise_accepted = true;
    if (aggregate_accepted && require_componentwise_capacity) {
        const Real capacity = density.target_density() * density.bin_area();
        for (const auto& [bin, old_occupancy] : saved_bins) {
            const Real old_excess = std::max<Real>(old_occupancy - capacity, 0.0);
            const Real new_excess = std::max<Real>(occupancy[bin] - capacity, 0.0);
            if (new_excess > old_excess + tolerance) {
                componentwise_accepted = false;
                break;
            }
        }
    }
    if (aggregate_accepted && componentwise_accepted && commit) {
        return AtomicMoveResult::Accepted;
    }
    for (std::size_t i = 0; i < members.size(); ++i) {
        db.nodes[members[i]].x = old_positions[i].first;
        db.nodes[members[i]].y = old_positions[i].second;
    }
    for (const auto& [bin, value] : saved_bins) occupancy[bin] = value;
    if (aggregate_accepted && componentwise_accepted) return AtomicMoveResult::Accepted;
    return aggregate_accepted ? AtomicMoveResult::CapacityRejected
                              : AtomicMoveResult::Rejected;
}

std::vector<std::pair<Real, Real>> net_targets(const Database& db) {
    std::vector<Real> sum_x(db.nodes.size(), 0.0);
    std::vector<Real> sum_y(db.nodes.size(), 0.0);
    std::vector<int> counts(db.nodes.size(), 0);
    for (const Net& net : db.nets) {
        if (net.pin_count == 0) continue;
        Real centroid_x = 0.0;
        Real centroid_y = 0.0;
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            const Pin& pin = db.pins[p];
            centroid_x += db.nodes[pin.node].x + pin.offset_x;
            centroid_y += db.nodes[pin.node].y + pin.offset_y;
        }
        centroid_x /= static_cast<Real>(net.pin_count);
        centroid_y /= static_cast<Real>(net.pin_count);
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            const int id = db.pins[p].node;
            sum_x[id] += centroid_x;
            sum_y[id] += centroid_y;
            ++counts[id];
        }
    }
    std::vector<std::pair<Real, Real>> result(db.nodes.size());
    for (std::size_t id = 0; id < db.nodes.size(); ++id) {
        if (counts[id] > 0) {
            result[id] = {sum_x[id] / counts[id], sum_y[id] / counts[id]};
        } else {
            result[id] = {db.nodes[id].x, db.nodes[id].y};
        }
    }
    return result;
}

struct Destination {
    Real distance = 0.0;
    int bin = -1;
};

struct AtomicCandidate {
    Real x = 0.0;
    Real y = 0.0;
    Real first_overflow_delta = 0.0;
    Real individual_score = 0.0;
    int bin = -1;
};

struct AtomicGroupBid {
    Real hpwl_delta = 0.0;
    int batch_group = -1;
    int group_id = -1;
    AtomicCandidate candidate;
};

struct FartherDestination {
    bool operator()(const Destination& a, const Destination& b) const noexcept {
        return a.distance > b.distance;
    }
};

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

std::vector<int> build_connectivity_rank(
    const Database& db, const std::vector<std::vector<int>>& node_nets,
    int degree_limit) {
    std::vector<int> rank(db.nodes.size(), std::numeric_limits<int>::max());
    std::vector<unsigned char> queued(db.nodes.size(), 0);
    std::vector<unsigned char> visited_net(db.nets.size(), 0);
    std::vector<int> stack;
    int next_rank = 0;
    for (int seed : db.movable_ids) {
        if (queued[seed]) continue;
        queued[seed] = 1;
        stack.push_back(seed);
        while (!stack.empty()) {
            const int id = stack.back();
            stack.pop_back();
            rank[id] = next_rank++;
            const auto& nets = node_nets[id];
            for (auto net_it = nets.rbegin(); net_it != nets.rend(); ++net_it) {
                const int net_id = *net_it;
                const Net& net = db.nets[net_id];
                if (visited_net[net_id] || net.pin_count >
                        static_cast<std::size_t>(degree_limit)) {
                    continue;
                }
                visited_net[net_id] = 1;
                for (std::size_t offset = net.pin_count; offset > 0; --offset) {
                    const int neighbor = db.pins[net.pin_begin + offset - 1].node;
                    if (!db.nodes[neighbor].fixed && !queued[neighbor]) {
                        queued[neighbor] = 1;
                        stack.push_back(neighbor);
                    }
                }
            }
        }
    }
    return rank;
}

Real incident_hpwl_delta(const Database& db, int node_id, Real new_x, Real new_y,
                         const std::vector<int>& nets, int degree_limit) {
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
            const Real old_x = node.x + pin.offset_x;
            const Real old_y = node.y + pin.offset_y;
            const Real candidate_x = (pin.node == node_id ? new_x : node.x) + pin.offset_x;
            const Real candidate_y = (pin.node == node_id ? new_y : node.y) + pin.offset_y;
            old_min_x = std::min(old_min_x, old_x);
            old_max_x = std::max(old_max_x, old_x);
            old_min_y = std::min(old_min_y, old_y);
            old_max_y = std::max(old_max_y, old_y);
            new_min_x = std::min(new_min_x, candidate_x);
            new_max_x = std::max(new_max_x, candidate_x);
            new_min_y = std::min(new_min_y, candidate_y);
            new_max_y = std::max(new_max_y, candidate_y);
        }
        delta += net.weight * ((new_max_x - new_min_x) + (new_max_y - new_min_y) -
                               (old_max_x - old_min_x) - (old_max_y - old_min_y));
    }
    return delta;
}

Real group_incident_hpwl_delta(
    const Database& db, const std::vector<int>& members, Real dx, Real dy,
    const std::vector<std::vector<int>>& node_nets) {
    std::vector<int> member_ids = members;
    std::sort(member_ids.begin(), member_ids.end());
    std::vector<int> nets;
    for (int id : members) {
        nets.insert(nets.end(), node_nets[id].begin(), node_nets[id].end());
    }
    std::sort(nets.begin(), nets.end());
    nets.erase(std::unique(nets.begin(), nets.end()), nets.end());

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
            const bool translated = std::binary_search(
                member_ids.begin(), member_ids.end(), pin.node);
            const Real new_x = old_x + (translated ? dx : 0.0);
            const Real new_y = old_y + (translated ? dy : 0.0);
            old_min_x = std::min(old_min_x, old_x);
            old_max_x = std::max(old_max_x, old_x);
            old_min_y = std::min(old_min_y, old_y);
            old_max_y = std::max(old_max_y, old_y);
            new_min_x = std::min(new_min_x, new_x);
            new_max_x = std::max(new_max_x, new_x);
            new_min_y = std::min(new_min_y, new_y);
            new_max_y = std::max(new_max_y, new_y);
        }
        result += net.weight * ((new_max_x - new_min_x) + (new_max_y - new_min_y) -
                                (old_max_x - old_min_x) - (old_max_y - old_min_y));
    }
    return result;
}

std::uint64_t splitmix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}

struct AuctionCandidate {
    Real x = 0.0;
    Real y = 0.0;
    Real target_distance = 0.0;
    MoveDelta density_delta;
};

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

struct HilbertMapping {
    std::vector<int> bins;
    std::vector<int> node_target_position;
};

HilbertMapping build_hilbert_mapping(
    const Database& db, const ExactOverlapDensity& density,
    const std::vector<int>& connectivity_rank) {
    if (density.bins_x() != density.bins_y() ||
        (density.bins_x() & (density.bins_x() - 1)) != 0) {
        throw std::invalid_argument("Hilbert transport requires a square power-of-two grid");
    }
    const int n = density.bins_x();
    std::vector<std::pair<std::uint64_t, int>> ordered;
    ordered.reserve(static_cast<std::size_t>(n) * n);
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            ordered.emplace_back(hilbert_index(n, x, y), y * n + x);
        }
    }
    std::sort(ordered.begin(), ordered.end());
    HilbertMapping result;
    result.bins.reserve(ordered.size());
    std::vector<Real> cumulative_capacity;
    cumulative_capacity.reserve(ordered.size());
    const Real bin_capacity = density.target_density() * density.bin_area();
    Real total_capacity = 0.0;
    for (const auto& entry : ordered) {
        const int bin = entry.second;
        result.bins.push_back(bin);
        total_capacity += std::max<Real>(
            bin_capacity - density.fixed_occupancy()[bin], 0.0);
        cumulative_capacity.push_back(total_capacity);
    }
    if (total_capacity <= 0.0 || db.movable_area <= 0.0) {
        throw std::runtime_error("Hilbert transport has no available capacity");
    }
    std::vector<int> nodes = db.movable_ids;
    std::sort(nodes.begin(), nodes.end(), [&](int a, int b) {
        return connectivity_rank[a] < connectivity_rank[b];
    });
    result.node_target_position.assign(db.nodes.size(), 0);
    Real cumulative_area = 0.0;
    for (int id : nodes) {
        const Real midpoint_area = cumulative_area + 0.5 * db.nodes[id].area();
        const Real target_capacity = midpoint_area / db.movable_area * total_capacity;
        const auto it = std::lower_bound(
            cumulative_capacity.begin(), cumulative_capacity.end(), target_capacity);
        result.node_target_position[id] = static_cast<int>(
            std::distance(cumulative_capacity.begin(), it));
        cumulative_area += db.nodes[id].area();
    }
    return result;
}

std::pair<std::vector<int>, int> build_transport_groups(
    const Database& db, int group_size, int degree_limit, bool weighted_grouping) {
    std::vector<int> group(db.nodes.size(), -1);
    int next_group = 0;
    if (group_size > 1 && weighted_grouping) {
        const auto node_nets = build_node_nets(db);
        std::vector<int> seeds = db.movable_ids;
        std::sort(seeds.begin(), seeds.end(), [&](int lhs, int rhs) {
            if (db.node_pin_count[lhs] != db.node_pin_count[rhs]) {
                return db.node_pin_count[lhs] > db.node_pin_count[rhs];
            }
            return lhs < rhs;
        });
        std::vector<Real> score(db.nodes.size(), 0.0);
        std::vector<int> stamp(db.nodes.size(), 0);
        int current_stamp = 0;
        for (int seed : seeds) {
            if (group[seed] >= 0) continue;
            ++current_stamp;
            std::vector<int> candidates;
            for (int net_id : node_nets[seed]) {
                const Net& net = db.nets[net_id];
                if (net.pin_count < 2 ||
                    net.pin_count > static_cast<std::size_t>(degree_limit)) {
                    continue;
                }
                std::vector<int> net_neighbors;
                net_neighbors.reserve(net.pin_count);
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const int neighbor = db.pins[p].node;
                    if (neighbor == seed || db.nodes[neighbor].fixed ||
                        group[neighbor] >= 0 ||
                        std::find(net_neighbors.begin(), net_neighbors.end(), neighbor) !=
                            net_neighbors.end()) {
                        continue;
                    }
                    net_neighbors.push_back(neighbor);
                }
                const Real contribution = net.weight /
                    static_cast<Real>(std::max<std::size_t>(net.pin_count - 1, 1));
                for (int neighbor : net_neighbors) {
                    if (stamp[neighbor] != current_stamp) {
                        stamp[neighbor] = current_stamp;
                        score[neighbor] = 0.0;
                        candidates.push_back(neighbor);
                    }
                    score[neighbor] += contribution;
                }
            }
            if (candidates.empty()) continue;
            std::sort(candidates.begin(), candidates.end(), [&](int lhs, int rhs) {
                if (score[lhs] != score[rhs]) return score[lhs] > score[rhs];
                return lhs < rhs;
            });
            group[seed] = next_group;
            const int take = std::min<int>(
                group_size - 1, static_cast<int>(candidates.size()));
            for (int i = 0; i < take; ++i) group[candidates[i]] = next_group;
            ++next_group;
        }
    } else if (group_size > 1) {
        std::vector<int> net_order;
        for (int net_id = 0; net_id < static_cast<int>(db.nets.size()); ++net_id) {
            const std::size_t degree = db.nets[net_id].pin_count;
            if (degree >= 2 && degree <= static_cast<std::size_t>(degree_limit)) {
                net_order.push_back(net_id);
            }
        }
        std::stable_sort(net_order.begin(), net_order.end(), [&](int a, int b) {
            return db.nets[a].pin_count < db.nets[b].pin_count;
        });
        for (int net_id : net_order) {
            const Net& net = db.nets[net_id];
            if (net.pin_count < 2 || net.pin_count > static_cast<std::size_t>(degree_limit)) {
                continue;
            }
            std::vector<int> candidates;
            candidates.reserve(net.pin_count);
            for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
                const int id = db.pins[p].node;
                if (!db.nodes[id].fixed && group[id] < 0 &&
                    std::find(candidates.begin(), candidates.end(), id) == candidates.end()) {
                    candidates.push_back(id);
                }
            }
            for (std::size_t begin = 0; begin + 1 < candidates.size();) {
                const std::size_t end = std::min(
                    candidates.size(), begin + static_cast<std::size_t>(group_size));
                if (end - begin < 2) break;
                for (std::size_t i = begin; i < end; ++i) group[candidates[i]] = next_group;
                ++next_group;
                begin = end;
            }
        }
    }
    for (int id : db.movable_ids) {
        if (group[id] < 0) group[id] = next_group++;
    }
    return {std::move(group), next_group};
}

struct IdentityBucketKey {
    std::int64_t width = 0;
    std::int64_t height = 0;
    int bin = 0;

    bool operator==(const IdentityBucketKey& other) const noexcept {
        return width == other.width && height == other.height && bin == other.bin;
    }
};

struct IdentityBucketHash {
    std::size_t operator()(const IdentityBucketKey& key) const noexcept {
        std::size_t value = std::hash<std::int64_t>{}(key.width);
        value ^= std::hash<std::int64_t>{}(key.height) + 0x9e3779b9U +
                 (value << 6U) + (value >> 2U);
        value ^= std::hash<int>{}(key.bin) + 0x9e3779b9U +
                 (value << 6U) + (value >> 2U);
        return value;
    }
};

std::int64_t identity_dimension_key(Real value) {
    return static_cast<std::int64_t>(std::llround(value * 1000000.0));
}

int identity_spatial_bin(const Database& db, const ExactOverlapDensity& density,
                         Real x, Real y) {
    const int bx = std::clamp(
        static_cast<int>(std::floor((x - db.xl) / density.bin_width())),
        0, density.bins_x() - 1);
    const int by = std::clamp(
        static_cast<int>(std::floor((y - db.yl) / density.bin_height())),
        0, density.bins_y() - 1);
    return by * density.bins_x() + bx;
}

Real exact_net_subset_hpwl(const Database& db, const std::vector<int>& nets) {
    Real result = 0.0;
    for (int net_id : nets) {
        const Net& net = db.nets[net_id];
        if (net.pin_count < 2) continue;
        Real min_x = std::numeric_limits<Real>::infinity();
        Real max_x = -std::numeric_limits<Real>::infinity();
        Real min_y = std::numeric_limits<Real>::infinity();
        Real max_y = -std::numeric_limits<Real>::infinity();
        for (std::size_t p = net.pin_begin; p < net.pin_begin + net.pin_count; ++p) {
            const Pin& pin = db.pins[p];
            const Node& node = db.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            min_x = std::min(min_x, x);
            max_x = std::max(max_x, x);
            min_y = std::min(min_y, y);
            max_y = std::max(max_y, y);
        }
        result += net.weight * ((max_x - min_x) + (max_y - min_y));
    }
    return result;
}

void run_identity_exchanges(
    Database& db, ExactOverlapDensity& density, const TransportConfig& config,
    const std::vector<int>& group, int group_count,
    const std::vector<std::vector<int>>& node_nets,
    const std::vector<Real>& occupancy, Real capacity, Real tolerance,
    TransportStats& stats) {
    if (config.identity_exchange_passes == 0) return;

    const auto started = std::chrono::steady_clock::now();
    ExactHpwl hpwl(db);
    stats.identity_initial_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    stats.identity_initial_overflow =
        density.evaluate(0.0, 1.0, nullptr, nullptr).overflow;
    Real audited_hpwl = stats.identity_initial_hpwl;
    Real audited_overflow = stats.identity_initial_overflow;

    for (int pass = 0; pass < config.identity_exchange_passes; ++pass) {
        const auto pass_targets = net_targets(db);
        std::unordered_map<IdentityBucketKey, std::vector<int>, IdentityBucketHash>
            buckets;
        buckets.reserve(db.movable_ids.size() * 2);
        for (int id : db.movable_ids) {
            const Node& node = db.nodes[id];
            buckets[{identity_dimension_key(node.width),
                     identity_dimension_key(node.height),
                     identity_spatial_bin(db, density, node.x, node.y)}]
                .push_back(id);
        }

        std::vector<std::vector<int>> contributors(occupancy.size());
        std::vector<AreaChange> footprint;
        for (int id : db.movable_ids) {
            footprint.clear();
            const Node& node = db.nodes[id];
            rectangle_changes(db, density, node, node.x, node.y, 1.0, footprint);
            for (const AreaChange& item : footprint) contributors[item.bin].push_back(id);
        }
        std::vector<int> sources;
        for (int bin = 0; bin < static_cast<int>(occupancy.size()); ++bin) {
            if (occupancy[bin] > capacity + tolerance && !contributors[bin].empty()) {
                sources.push_back(bin);
            }
        }
        std::sort(sources.begin(), sources.end(), [&](int a, int b) {
            const Real excess_a = occupancy[a] - capacity;
            const Real excess_b = occupancy[b] - capacity;
            return excess_a == excess_b ? a < b : excess_a > excess_b;
        });
        if (static_cast<int>(sources.size()) > config.max_source_bins) {
            sources.resize(config.max_source_bins);
        }

        std::vector<unsigned char> used(db.nodes.size(), 0);
        std::vector<Real> pass_old_x(db.movable_ids.size());
        std::vector<Real> pass_old_y(db.movable_ids.size());
        for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
            const Node& node = db.nodes[db.movable_ids[i]];
            pass_old_x[i] = node.x;
            pass_old_y[i] = node.y;
        }

        for (int source : sources) {
            std::vector<std::vector<int>> active(group_count);
            for (int id : contributors[source]) {
                if (!used[id]) active[group[id]].push_back(id);
            }
            const int source_x = source % density.bins_x();
            const int source_y = source / density.bins_x();

            for (int group_id = 0; group_id < group_count; ++group_id) {
                std::vector<int>& members = active[group_id];
                members.erase(std::remove_if(members.begin(), members.end(),
                                             [&](int id) { return used[id] != 0; }),
                              members.end());
                std::sort(members.begin(), members.end());
                members.erase(std::unique(members.begin(), members.end()), members.end());
                if (members.size() < 2) continue;

                Real mean_target_x = 0.0;
                Real mean_target_y = 0.0;
                for (int id : members) {
                    mean_target_x += pass_targets[id].first;
                    mean_target_y += pass_targets[id].second;
                }
                mean_target_x /= static_cast<Real>(members.size());
                mean_target_y /= static_cast<Real>(members.size());

                std::vector<int> destination_bins;
                const auto append_destination = [&](Real x, Real y) {
                    const int bin = identity_spatial_bin(db, density, x, y);
                    if (std::find(destination_bins.begin(), destination_bins.end(), bin) ==
                        destination_bins.end()) {
                        destination_bins.push_back(bin);
                    }
                };
                append_destination(mean_target_x, mean_target_y);
                std::vector<int> target_order = members;
                std::sort(target_order.begin(), target_order.end(), [&](int a, int b) {
                    const Real adx = pass_targets[a].first - mean_target_x;
                    const Real ady = pass_targets[a].second - mean_target_y;
                    const Real bdx = pass_targets[b].first - mean_target_x;
                    const Real bdy = pass_targets[b].second - mean_target_y;
                    const Real da = adx * adx + ady * ady;
                    const Real db_value = bdx * bdx + bdy * bdy;
                    return da == db_value ? a < b : da < db_value;
                });
                for (int id : target_order) append_destination(
                    pass_targets[id].first, pass_targets[id].second);
                if (static_cast<int>(destination_bins.size()) >
                    config.identity_exchange_candidates) {
                    destination_bins.resize(config.identity_exchange_candidates);
                }

                Real best_delta = 0.0;
                std::vector<int> best_donors;
                for (int destination : destination_bins) {
                    const int destination_x = destination % density.bins_x();
                    const int destination_y = destination / density.bins_x();
                    const Real dx = (destination_x - source_x) * density.bin_width();
                    const Real dy = (destination_y - source_y) * density.bin_height();
                    std::vector<int> donors;
                    donors.reserve(members.size());
                    bool complete = true;
                    for (int member : members) {
                        const Node& member_node = db.nodes[member];
                        const Real translated_x = member_node.x + dx;
                        const Real translated_y = member_node.y + dy;
                        const int target_bin = identity_spatial_bin(
                            db, density, translated_x, translated_y);
                        const int target_x = target_bin % density.bins_x();
                        const int target_y = target_bin / density.bins_x();
                        Real best_distance = std::numeric_limits<Real>::infinity();
                        int best_donor = -1;
                        for (int radius = 0; radius <= config.identity_exchange_radius;
                             ++radius) {
                            const int x0 = std::max(0, target_x - radius);
                            const int x1 = std::min(density.bins_x() - 1,
                                                    target_x + radius);
                            const int y0 = std::max(0, target_y - radius);
                            const int y1 = std::min(density.bins_y() - 1,
                                                    target_y + radius);
                            for (int by = y0; by <= y1; ++by) {
                                for (int bx = x0; bx <= x1; ++bx) {
                                    if (radius > 0 && bx > x0 && bx < x1 &&
                                        by > y0 && by < y1) {
                                        continue;
                                    }
                                    const IdentityBucketKey key{
                                        identity_dimension_key(member_node.width),
                                        identity_dimension_key(member_node.height),
                                        by * density.bins_x() + bx};
                                    const auto found = buckets.find(key);
                                    if (found == buckets.end()) continue;
                                    for (int donor : found->second) {
                                        if (used[donor] || donor == member ||
                                            std::find(members.begin(), members.end(), donor) !=
                                                members.end() ||
                                            std::find(donors.begin(), donors.end(), donor) !=
                                                donors.end()) {
                                            continue;
                                        }
                                        const Node& donor_node = db.nodes[donor];
                                        if (donor_node.fixed ||
                                            donor_node.width != member_node.width ||
                                            donor_node.height != member_node.height ||
                                            identity_spatial_bin(db, density, donor_node.x,
                                                                 donor_node.y) != key.bin) {
                                            continue;
                                        }
                                        const Real ddx = donor_node.x - translated_x;
                                        const Real ddy = donor_node.y - translated_y;
                                        const Real distance = ddx * ddx + ddy * ddy;
                                        if (distance < best_distance ||
                                            (distance == best_distance && donor < best_donor)) {
                                            best_distance = distance;
                                            best_donor = donor;
                                        }
                                    }
                                }
                            }
                        }
                        if (best_donor < 0) {
                            complete = false;
                            break;
                        }
                        donors.push_back(best_donor);
                    }
                    if (!complete) continue;

                    std::vector<int> nets;
                    for (int id : members) {
                        nets.insert(nets.end(), node_nets[id].begin(), node_nets[id].end());
                    }
                    for (int id : donors) {
                        nets.insert(nets.end(), node_nets[id].begin(), node_nets[id].end());
                    }
                    std::sort(nets.begin(), nets.end());
                    nets.erase(std::unique(nets.begin(), nets.end()), nets.end());
                    const Real before = exact_net_subset_hpwl(db, nets);
                    for (std::size_t i = 0; i < members.size(); ++i) {
                        std::swap(db.nodes[members[i]].x, db.nodes[donors[i]].x);
                        std::swap(db.nodes[members[i]].y, db.nodes[donors[i]].y);
                    }
                    const Real after = exact_net_subset_hpwl(db, nets);
                    for (std::size_t i = 0; i < members.size(); ++i) {
                        std::swap(db.nodes[members[i]].x, db.nodes[donors[i]].x);
                        std::swap(db.nodes[members[i]].y, db.nodes[donors[i]].y);
                    }
                    ++stats.identity_exchange_attempts;
                    const Real delta = after - before;
                    if (delta < best_delta - 1.0e-12) {
                        best_delta = delta;
                        best_donors = std::move(donors);
                    }
                }
                if (best_donors.empty()) continue;
                for (std::size_t i = 0; i < members.size(); ++i) {
                    std::swap(db.nodes[members[i]].x, db.nodes[best_donors[i]].x);
                    std::swap(db.nodes[members[i]].y, db.nodes[best_donors[i]].y);
                    used[members[i]] = 1;
                    used[best_donors[i]] = 1;
                }
                ++stats.identity_exchange_accepted;
                stats.identity_exchange_permuted_nodes +=
                    2 * static_cast<int>(members.size());
            }
        }

        const Real pass_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        const Real pass_overflow =
            density.evaluate(0.0, 1.0, nullptr, nullptr).overflow;
        if (pass_hpwl > audited_hpwl + 1.0e-9 ||
            std::abs(pass_overflow - audited_overflow) > 1.0e-12) {
            for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
                Node& node = db.nodes[db.movable_ids[i]];
                node.x = pass_old_x[i];
                node.y = pass_old_y[i];
            }
            throw std::logic_error("identity exchange failed its exact-oracle audit");
        }
        audited_hpwl = pass_hpwl;
        audited_overflow = pass_overflow;
    }

    stats.identity_final_hpwl = audited_hpwl;
    stats.identity_final_overflow = audited_overflow;
    stats.identity_wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
}

}  // namespace

TransportDestinationMode parse_transport_destination_mode(const std::string& name) {
    if (name == "nearest") return TransportDestinationMode::Nearest;
    if (name == "auction") return TransportDestinationMode::Auction;
    if (name == "hilbert") return TransportDestinationMode::Hilbert;
    throw std::invalid_argument("unknown transport destination mode: " + name);
}

const char* transport_destination_mode_name(TransportDestinationMode mode) noexcept {
    if (mode == TransportDestinationMode::Auction) return "auction";
    if (mode == TransportDestinationMode::Hilbert) return "hilbert";
    return "nearest";
}

TransportStats transport_excess_to_capacity(
    Database& db, ExactOverlapDensity& density, const TransportConfig& config) {
    if (config.rounds < 0 || config.max_moves < 0 || config.max_source_bins <= 0 ||
        config.candidate_lookahead <= 0 || config.group_size <= 0 ||
        config.group_degree_limit < 2 || config.group_destination_radius < 0 ||
        config.identity_exchange_passes < 0 ||
        config.identity_exchange_candidates <= 0 ||
        config.identity_exchange_radius < 0 ||
        config.rigid_rounds < -1 || config.atomic_rounds < -1 ||
        config.atomic_capacity_candidates < 1 ||
        (config.atomic_capacity_candidates > 1 &&
         !config.atomic_componentwise_capacity) ||
        (config.atomic_source_batch &&
         (!config.atomic_source_active || config.atomic_capacity_candidates < 2)) ||
        config.atomic_min_gain_ratio < 0.0 || config.atomic_min_gain_ratio > 1.0 ||
        config.auction_candidates <= 0 ||
        config.auction_shortlist <= 0 ||
        config.auction_shortlist > config.auction_candidates ||
        config.auction_hpwl_degree_limit < 2 ||
        config.auction_density_price_weight < 0.0 || config.auction_rounds < -1 ||
        config.connectivity_degree_limit < 2 || config.hilbert_rounds < 1 ||
        config.hilbert_window < 1 || config.global_hpwl_budget < -1.0) {
        throw std::invalid_argument("invalid capacity transport configuration");
    }
    TransportStats stats;
    if (config.rounds == 0 && config.identity_exchange_passes == 0) return stats;

    const auto started = std::chrono::steady_clock::now();
    ExactHpwl hpwl(db);
    const DensityMetrics initial_density = density.evaluate(0.0, 1.0, nullptr, nullptr);
    stats.initial_overflow = initial_density.overflow;
    stats.initial_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    std::vector<Real> occupancy = density.occupancy();
    const Real capacity = density.target_density() * density.bin_area();
    const Real tolerance = 1.0e-10 * density.bin_area();
    std::vector<unsigned char> moved(db.nodes.size(), 0);
    auto [group, group_count] = build_transport_groups(
        db, config.group_size, config.group_degree_limit, config.weighted_grouping);
    stats.groups = group_count;
    std::vector<int> group_members(group_count, 0);
    for (int id : db.movable_ids) ++group_members[group[id]];
    std::vector<std::vector<int>> group_nodes(group_count);
    for (int id : db.movable_ids) group_nodes[group[id]].push_back(id);
    std::vector<int> group_destination(group_count, -1);
    std::vector<std::pair<Real, Real>> group_targets(group_count, {0.0, 0.0});
    std::vector<Real> group_dx(group_count, 0.0);
    std::vector<Real> group_dy(group_count, 0.0);
    std::vector<unsigned char> group_has_displacement(group_count, 0);
    std::vector<unsigned char> atomic_group_processed(group_count, 0);
    std::vector<unsigned char> atomic_group_fallback(group_count, 0);
    std::vector<int> atomic_source_stamp(group_count, 0);
    int source_stamp = 0;
    const auto node_nets = build_node_nets(db);
    const auto connectivity_rank =
        config.connectivity_order ||
                config.destination_mode == TransportDestinationMode::Hilbert
        ? build_connectivity_rank(db, node_nets, config.connectivity_degree_limit)
        : std::vector<int>(db.nodes.size(), 0);
    const HilbertMapping hilbert =
        config.destination_mode == TransportDestinationMode::Hilbert
            ? build_hilbert_mapping(db, density, connectivity_rank)
            : HilbertMapping{};

    run_identity_exchanges(db, density, config, group, group_count, node_nets,
                           occupancy, capacity, tolerance, stats);

    for (int round = 0; round < config.rounds && stats.moves < config.max_moves; ++round) {
        const Real round_start_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        std::vector<std::pair<int, std::pair<Real, Real>>> round_positions;
        if (config.global_hpwl_budget >= 0.0) {
            round_positions.reserve(db.movable_ids.size());
            for (int id : db.movable_ids) {
                round_positions.push_back({id, {db.nodes[id].x, db.nodes[id].y}});
            }
        }
        const std::vector<Real> round_occupancy = occupancy;
        const int moves_before_round = stats.moves;
        TransportDestinationMode round_mode = config.destination_mode;
        if (round_mode == TransportDestinationMode::Auction &&
            config.auction_rounds >= 0 && round >= config.auction_rounds) {
            round_mode = TransportDestinationMode::Nearest;
        }
        if (round_mode == TransportDestinationMode::Hilbert &&
            round >= config.hilbert_rounds) {
            round_mode = TransportDestinationMode::Nearest;
        }
        const bool rigid_groups = config.rigid_group_displacement &&
            (config.rigid_rounds < 0 || round < config.rigid_rounds);
        const bool atomic_enabled = config.atomic_groups &&
            (config.atomic_rounds < 0 || round < config.atomic_rounds);
        std::fill(moved.begin(), moved.end(), 0);
        std::fill(group_has_displacement.begin(), group_has_displacement.end(), 0);
        std::fill(atomic_group_processed.begin(), atomic_group_processed.end(), 0);
        std::fill(atomic_group_fallback.begin(), atomic_group_fallback.end(), 0);
        const auto targets = net_targets(db);
        if (config.group_collective_target) {
            std::fill(group_targets.begin(), group_targets.end(),
                      std::pair<Real, Real>{0.0, 0.0});
            for (int id : db.movable_ids) {
                group_targets[group[id]].first += targets[id].first;
                group_targets[group[id]].second += targets[id].second;
            }
            for (int group_id = 0; group_id < group_count; ++group_id) {
                group_targets[group_id].first /= group_members[group_id];
                group_targets[group_id].second /= group_members[group_id];
            }
        }
        std::vector<std::vector<int>> contributors(occupancy.size());
        std::vector<AreaChange> footprint;
        for (int id : db.movable_ids) {
            footprint.clear();
            const Node& node = db.nodes[id];
            rectangle_changes(db, density, node, node.x, node.y, 1.0, footprint);
            for (const AreaChange& item : footprint) contributors[item.bin].push_back(id);
        }

        std::vector<int> sources;
        sources.reserve(occupancy.size());
        for (int bin = 0; bin < static_cast<int>(occupancy.size()); ++bin) {
            if (occupancy[bin] > capacity + tolerance && !contributors[bin].empty()) {
                sources.push_back(bin);
            }
        }
        std::sort(sources.begin(), sources.end(), [&](int a, int b) {
            return occupancy[a] - capacity > occupancy[b] - capacity;
        });
        if (static_cast<int>(sources.size()) > config.max_source_bins) {
            sources.resize(config.max_source_bins);
        }

        int round_moves = 0;
        for (int source : sources) {
            if (stats.moves >= config.max_moves) break;
            ++source_stamp;
            auto& nodes = contributors[source];
            std::sort(nodes.begin(), nodes.end(), [&](int a, int b) {
                if (config.connectivity_order && connectivity_rank[a] != connectivity_rank[b]) {
                    return connectivity_rank[a] < connectivity_rank[b];
                }
                if (group[a] != group[b]) return group[a] < group[b];
                if (db.node_pin_count[a] != db.node_pin_count[b]) {
                    return db.node_pin_count[a] < db.node_pin_count[b];
                }
                return db.nodes[a].area() < db.nodes[b].area();
            });
            const int source_y = source / density.bins_x();
            const int source_x = source % density.bins_x();
            std::priority_queue<Destination, std::vector<Destination>, FartherDestination> heap;
            if (round_mode == TransportDestinationMode::Nearest) {
                for (int by = 0; by < density.bins_y(); ++by) {
                    for (int bx = 0; bx < density.bins_x(); ++bx) {
                        const int bin = by * density.bins_x() + bx;
                        if (occupancy[bin] >= capacity - tolerance) continue;
                        const Real dx = static_cast<Real>(bx - source_x);
                        const Real dy = static_cast<Real>(by - source_y);
                        heap.push({dx * dx + dy * dy, bin});
                    }
                }
            }

            const auto coordinates_for = [&](const Node& candidate_node, int bx, int by) {
                Real x = db.xl + (bx + 0.5) * density.bin_width();
                Real y = db.yl + (by + 0.5) * density.bin_height();
                if (config.preserve_bin_offset) {
                    x = candidate_node.x + (bx - source_x) * density.bin_width();
                    y = candidate_node.y + (by - source_y) * density.bin_height();
                }
                x = std::clamp(x, db.xl + 0.5 * candidate_node.width,
                               db.xh - 0.5 * candidate_node.width);
                y = std::clamp(y, db.yl + 0.5 * candidate_node.height,
                               db.yh - 0.5 * candidate_node.height);
                return std::pair<Real, Real>{x, y};
            };

            const bool source_batch = round_mode == TransportDestinationMode::Nearest &&
                atomic_enabled && config.atomic_source_batch;
            if (source_batch) {
                std::vector<int> batch_group_ids;
                std::vector<std::vector<int>> batch_members;
                for (int id : nodes) {
                    if (moved[id]) continue;
                    const int group_id = group[id];
                    if (group_members[group_id] < 2 || atomic_group_fallback[group_id]) {
                        continue;
                    }
                    if (batch_group_ids.empty() || batch_group_ids.back() != group_id) {
                        batch_group_ids.push_back(group_id);
                        batch_members.emplace_back();
                    }
                    batch_members.back().push_back(id);
                }

                std::vector<Destination> global_destinations;
                auto batch_heap = heap;
                while (!batch_heap.empty() && static_cast<int>(global_destinations.size()) <
                                                 config.candidate_lookahead) {
                    global_destinations.push_back(batch_heap.top());
                    batch_heap.pop();
                }
                std::vector<AtomicGroupBid> bids;
                std::vector<unsigned char> eligible(batch_members.size(), 0);
                for (int batch_id = 0; batch_id < static_cast<int>(batch_members.size());
                     ++batch_id) {
                    auto& members = batch_members[batch_id];
                    if (members.size() < 2) continue;
                    eligible[batch_id] = 1;
                    ++stats.atomic_attempts;
                    ++stats.atomic_batch_groups;
                    const int group_id = batch_group_ids[batch_id];
                    const int first_id = members.front();
                    const Node& first = db.nodes[first_id];
                    const auto& score_target = config.group_collective_target
                        ? group_targets[group_id] : targets[first_id];

                    std::vector<Destination> destinations;
                    bool used_local_destinations = false;
                    if (config.group_destination_radius > 0 &&
                        group_destination[group_id] >= 0) {
                        used_local_destinations = true;
                        const int anchor = group_destination[group_id];
                        const int anchor_y = anchor / density.bins_x();
                        const int anchor_x = anchor % density.bins_x();
                        const int radius = config.group_destination_radius;
                        for (int by = std::max(0, anchor_y - radius);
                             by <= std::min(density.bins_y() - 1, anchor_y + radius); ++by) {
                            for (int bx = std::max(0, anchor_x - radius);
                                 bx <= std::min(density.bins_x() - 1, anchor_x + radius); ++bx) {
                                const int bin = by * density.bins_x() + bx;
                                if (occupancy[bin] >= capacity - tolerance) continue;
                                const Real dx = static_cast<Real>(bx - anchor_x);
                                const Real dy = static_cast<Real>(by - anchor_y);
                                destinations.push_back({dx * dx + dy * dy, bin});
                            }
                        }
                        std::sort(destinations.begin(), destinations.end(),
                                  [](const Destination& a, const Destination& b) {
                            if (a.distance != b.distance) return a.distance < b.distance;
                            return a.bin < b.bin;
                        });
                        if (static_cast<int>(destinations.size()) >
                            config.candidate_lookahead) {
                            destinations.resize(config.candidate_lookahead);
                        }
                    } else {
                        destinations = global_destinations;
                    }

                    std::vector<AtomicCandidate> candidates;
                    const auto append_candidates = [&](const std::vector<Destination>& pool) {
                        for (const Destination& destination : pool) {
                            const int by = destination.bin / density.bins_x();
                            const int bx = destination.bin % density.bins_x();
                            const auto [x, y] = coordinates_for(first, bx, by);
                            MoveDelta delta = evaluate_move(
                                db, density, first, x, y, occupancy);
                            if (delta.overflow_area >= -tolerance) continue;
                            const Real tx = x - score_target.first;
                            const Real ty = y - score_target.second;
                            const Real score = delta.energy +
                                1.0e-9 * (tx * tx + ty * ty);
                            candidates.push_back(
                                {x, y, delta.overflow_area, score, destination.bin});
                        }
                    };
                    append_candidates(destinations);
                    if (candidates.empty() && used_local_destinations) {
                        append_candidates(global_destinations);
                    }
                    std::sort(candidates.begin(), candidates.end(),
                              [](const AtomicCandidate& a, const AtomicCandidate& b) {
                        if (a.individual_score != b.individual_score) {
                            return a.individual_score < b.individual_score;
                        }
                        return a.bin < b.bin;
                    });
                    if (static_cast<int>(candidates.size()) >
                        config.atomic_capacity_candidates) {
                        candidates.resize(config.atomic_capacity_candidates);
                    }
                    for (const AtomicCandidate& candidate : candidates) {
                        const Real dx = candidate.x - first.x;
                        const Real dy = candidate.y - first.y;
                        const Real maximum_group_delta = config.atomic_min_gain_ratio *
                            candidate.first_overflow_delta;
                        const AtomicMoveResult trial = apply_atomic_group_move(
                            db, density, members, dx, dy, maximum_group_delta,
                            tolerance, true, false, occupancy);
                        ++stats.atomic_bid_trials;
                        if (trial == AtomicMoveResult::CapacityRejected) {
                            ++stats.atomic_capacity_rejections;
                        }
                        if (trial != AtomicMoveResult::Accepted) continue;
                        ++stats.atomic_feasible_bids;
                        const Real hpwl_delta = group_incident_hpwl_delta(
                            db, members, dx, dy, node_nets);
                        bids.push_back({hpwl_delta, batch_id, group_id, candidate});
                    }
                }

                if (!batch_group_ids.empty()) ++stats.atomic_batch_sources;
                stats.atomic_batch_bids += static_cast<int>(bids.size());
                std::sort(bids.begin(), bids.end(), [](const AtomicGroupBid& a,
                                                       const AtomicGroupBid& b) {
                    if (a.hpwl_delta != b.hpwl_delta) return a.hpwl_delta < b.hpwl_delta;
                    if (a.group_id != b.group_id) return a.group_id < b.group_id;
                    return a.candidate.bin < b.candidate.bin;
                });
                std::vector<unsigned char> assigned(batch_members.size(), 0);
                for (const AtomicGroupBid& bid : bids) {
                    if (assigned[bid.batch_group]) continue;
                    auto& members = batch_members[bid.batch_group];
                    if (stats.moves + static_cast<int>(members.size()) > config.max_moves) {
                        continue;
                    }
                    const Node& first = db.nodes[members.front()];
                    const Real dx = bid.candidate.x - first.x;
                    const Real dy = bid.candidate.y - first.y;
                    const Real maximum_group_delta = config.atomic_min_gain_ratio *
                        bid.candidate.first_overflow_delta;
                    const AtomicMoveResult committed = apply_atomic_group_move(
                        db, density, members, dx, dy, maximum_group_delta,
                        tolerance, true, true, occupancy);
                    if (committed != AtomicMoveResult::Accepted) {
                        ++stats.atomic_batch_rejections;
                        continue;
                    }
                    assigned[bid.batch_group] = 1;
                    ++stats.atomic_groups;
                    ++stats.atomic_batch_commits;
                    const int group_id = bid.group_id;
                    if (config.group_destination_radius > 0 &&
                        group_destination[group_id] < 0) {
                        group_destination[group_id] = bid.candidate.bin;
                    }
                    for (int member : members) moved[member] = 1;
                    const int accepted_nodes = static_cast<int>(members.size());
                    stats.moves += accepted_nodes;
                    round_moves += accepted_nodes;
                }
                for (int batch_id = 0; batch_id < static_cast<int>(batch_members.size());
                     ++batch_id) {
                    if (!eligible[batch_id] || assigned[batch_id]) continue;
                    atomic_group_fallback[batch_group_ids[batch_id]] = 1;
                    ++stats.atomic_fallbacks;
                }
            }

            for (int id : nodes) {
                if (stats.moves >= config.max_moves || occupancy[source] <= capacity + tolerance) break;
                if (moved[id]) continue;
                Node& node = db.nodes[id];
                const int group_id = group[id];
                std::vector<int> source_active_members;
                const std::vector<int>* atomic_members = &group_nodes[group_id];
                bool atomic_group = round_mode == TransportDestinationMode::Nearest &&
                    atomic_enabled && group_members[group_id] > 1 &&
                    !atomic_group_fallback[group_id];
                if (atomic_group && config.atomic_source_active) {
                    source_active_members.reserve(group_members[group_id]);
                    for (int candidate : nodes) {
                        if (!moved[candidate] && group[candidate] == group_id) {
                            source_active_members.push_back(candidate);
                        }
                    }
                    atomic_group = source_active_members.size() > 1;
                    atomic_members = &source_active_members;
                }
                if (atomic_group) {
                    if (config.atomic_source_active) {
                        if (atomic_source_stamp[group_id] == source_stamp) continue;
                        atomic_source_stamp[group_id] = source_stamp;
                    } else {
                        if (atomic_group_processed[group_id]) continue;
                        atomic_group_processed[group_id] = 1;
                    }
                }
                if (round_mode == TransportDestinationMode::Nearest &&
                    !atomic_group && rigid_groups && group_members[group_id] > 1 &&
                    group_has_displacement[group_id]) {
                    const Real x = node.x + group_dx[group_id];
                    const Real y = node.y + group_dy[group_id];
                    if (x < db.xl + 0.5 * node.width ||
                        x > db.xh - 0.5 * node.width ||
                        y < db.yl + 0.5 * node.height ||
                        y > db.yh - 0.5 * node.height) {
                        continue;
                    }
                    MoveDelta delta = evaluate_move(db, density, node, x, y, occupancy);
                    if (delta.overflow_area >= -tolerance) continue;
                    apply_move(node, x, y, delta, occupancy);
                    moved[id] = 1;
                    ++stats.moves;
                    ++round_moves;
                    continue;
                }
                if (round_mode == TransportDestinationMode::Hilbert) {
                    const int target = hilbert.node_target_position[id];
                    bool found = false;
                    Real best_x = node.x;
                    Real best_y = node.y;
                    MoveDelta best_delta;
                    for (int attempt = 0; attempt < config.hilbert_window; ++attempt) {
                        const int magnitude = (attempt + 1) / 2;
                        const int signed_offset = attempt == 0
                            ? 0
                            : (attempt % 2 == 1 ? magnitude : -magnitude);
                        const int position = target + signed_offset;
                        if (position < 0 || position >= static_cast<int>(hilbert.bins.size())) {
                            continue;
                        }
                        const int destination = hilbert.bins[position];
                        const int by = destination / density.bins_x();
                        const int bx = destination % density.bins_x();
                        const Real x = std::clamp(
                            db.xl + (bx + 0.5) * density.bin_width(),
                            db.xl + 0.5 * node.width, db.xh - 0.5 * node.width);
                        const Real y = std::clamp(
                            db.yl + (by + 0.5) * density.bin_height(),
                            db.yl + 0.5 * node.height, db.yh - 0.5 * node.height);
                        MoveDelta delta = evaluate_move(db, density, node, x, y, occupancy);
                        if (delta.overflow_area >= -tolerance) continue;
                        found = true;
                        best_x = x;
                        best_y = y;
                        best_delta = std::move(delta);
                        break;
                    }
                    if (!found) continue;
                    apply_move(node, best_x, best_y, best_delta, occupancy);
                    moved[id] = 1;
                    ++stats.moves;
                    ++round_moves;
                    continue;
                }
                if (round_mode == TransportDestinationMode::Auction) {
                    std::vector<AuctionCandidate> shortlist;
                    shortlist.reserve(config.auction_shortlist);
                    const std::uint64_t base = splitmix64(
                        static_cast<std::uint64_t>(group[id]) ^
                        (static_cast<std::uint64_t>(round + 1) << 32U));
                    for (int attempt = 0; attempt < config.auction_candidates; ++attempt) {
                        const std::uint64_t random = splitmix64(
                            base + static_cast<std::uint64_t>(attempt));
                        const int destination = static_cast<int>(
                            random % static_cast<std::uint64_t>(occupancy.size()));
                        const int by = destination / density.bins_x();
                        const int bx = destination % density.bins_x();
                        const Real x = std::clamp(
                            db.xl + (bx + 0.5) * density.bin_width(),
                            db.xl + 0.5 * node.width, db.xh - 0.5 * node.width);
                        const Real y = std::clamp(
                            db.yl + (by + 0.5) * density.bin_height(),
                            db.yl + 0.5 * node.height, db.yh - 0.5 * node.height);
                        MoveDelta delta = evaluate_move(db, density, node, x, y, occupancy);
                        if (delta.overflow_area >= -tolerance) continue;
                        const Real dx = x - targets[id].first;
                        const Real dy = y - targets[id].second;
                        AuctionCandidate candidate{x, y, dx * dx + dy * dy,
                                                   std::move(delta)};
                        if (static_cast<int>(shortlist.size()) < config.auction_shortlist) {
                            shortlist.push_back(std::move(candidate));
                        } else {
                            auto worst = std::max_element(
                                shortlist.begin(), shortlist.end(),
                                [](const AuctionCandidate& a, const AuctionCandidate& b) {
                                    return a.target_distance < b.target_distance;
                                });
                            if (candidate.target_distance < worst->target_distance) {
                                *worst = std::move(candidate);
                            }
                        }
                    }
                    bool found = false;
                    Real best_score = std::numeric_limits<Real>::infinity();
                    Real best_x = node.x;
                    Real best_y = node.y;
                    MoveDelta best_delta;
                    for (AuctionCandidate& candidate : shortlist) {
                        const Real hpwl_delta = incident_hpwl_delta(
                            db, id, candidate.x, candidate.y, node_nets[id],
                            config.auction_hpwl_degree_limit);
                        const Real score = hpwl_delta +
                            config.auction_density_price_weight * candidate.density_delta.energy;
                        if (!found || score < best_score) {
                            found = true;
                            best_score = score;
                            best_x = candidate.x;
                            best_y = candidate.y;
                            best_delta = std::move(candidate.density_delta);
                        }
                    }
                    if (!found) continue;
                    apply_move(node, best_x, best_y, best_delta, occupancy);
                    moved[id] = 1;
                    ++stats.moves;
                    ++round_moves;
                    continue;
                }
                std::vector<Destination> reusable;
                reusable.reserve(config.candidate_lookahead);
                bool found = false;
                Real best_x = node.x;
                Real best_y = node.y;
                Real best_score = std::numeric_limits<Real>::infinity();
                int best_bin = -1;
                MoveDelta best_delta;
                std::vector<AtomicCandidate> atomic_candidates;
                const auto& score_target =
                    config.group_collective_target && group_members[group_id] > 1
                    ? group_targets[group_id]
                    : targets[id];
                const auto candidate_coordinates = [&](int bx, int by) {
                    return coordinates_for(node, bx, by);
                };
                if (config.group_destination_radius > 0 &&
                    group_members[group_id] > 1 && group_destination[group_id] >= 0) {
                    const int anchor = group_destination[group_id];
                    const int anchor_y = anchor / density.bins_x();
                    const int anchor_x = anchor % density.bins_x();
                    std::vector<Destination> local;
                    const int radius = config.group_destination_radius;
                    for (int by = std::max(0, anchor_y - radius);
                         by <= std::min(density.bins_y() - 1, anchor_y + radius); ++by) {
                        for (int bx = std::max(0, anchor_x - radius);
                             bx <= std::min(density.bins_x() - 1, anchor_x + radius); ++bx) {
                            const int bin = by * density.bins_x() + bx;
                            if (occupancy[bin] >= capacity - tolerance) continue;
                            const Real dx = static_cast<Real>(bx - anchor_x);
                            const Real dy = static_cast<Real>(by - anchor_y);
                            local.push_back({dx * dx + dy * dy, bin});
                        }
                    }
                    std::sort(local.begin(), local.end(), [](const Destination& a,
                                                             const Destination& b) {
                        if (a.distance != b.distance) return a.distance < b.distance;
                        return a.bin < b.bin;
                    });
                    if (static_cast<int>(local.size()) > config.candidate_lookahead) {
                        local.resize(config.candidate_lookahead);
                    }
                    for (const Destination& destination : local) {
                        const int by = destination.bin / density.bins_x();
                        const int bx = destination.bin % density.bins_x();
                        const auto [x, y] = candidate_coordinates(bx, by);
                        MoveDelta delta = evaluate_move(db, density, node, x, y, occupancy);
                        if (delta.overflow_area >= -tolerance) continue;
                        const Real tx = x - score_target.first;
                        const Real ty = y - score_target.second;
                        const Real score = delta.energy + 1.0e-9 * (tx * tx + ty * ty);
                        if (atomic_group && config.atomic_capacity_candidates > 1) {
                            atomic_candidates.push_back(
                                {x, y, delta.overflow_area, score, destination.bin});
                        }
                        if (!found || score < best_score) {
                            found = true;
                            best_score = score;
                            best_x = x;
                            best_y = y;
                            best_bin = destination.bin;
                            best_delta = std::move(delta);
                        }
                    }
                }
                if (!found) {
                    for (int attempt = 0; attempt < config.candidate_lookahead &&
                                          !heap.empty(); ++attempt) {
                        const Destination destination = heap.top();
                        heap.pop();
                        const int by = destination.bin / density.bins_x();
                        const int bx = destination.bin % density.bins_x();
                        const auto [x, y] = candidate_coordinates(bx, by);
                        MoveDelta delta = evaluate_move(db, density, node, x, y, occupancy);
                        if (delta.overflow_area >= -tolerance) continue;
                        reusable.push_back(destination);
                        const Real tx = x - score_target.first;
                        const Real ty = y - score_target.second;
                        const Real target_distance = tx * tx + ty * ty;
                        const Real score = delta.energy + 1.0e-9 * target_distance;
                        if (atomic_group && config.atomic_capacity_candidates > 1) {
                            atomic_candidates.push_back(
                                {x, y, delta.overflow_area, score, destination.bin});
                        }
                        if (!found || score < best_score) {
                            found = true;
                            best_score = score;
                            best_x = x;
                            best_y = y;
                            best_bin = destination.bin;
                            best_delta = std::move(delta);
                        }
                    }
                }
                if (!found) {
                    if (atomic_group && config.atomic_fallback_individual) {
                        atomic_group_fallback[group_id] = 1;
                        ++stats.atomic_fallbacks;
                    }
                    continue;
                }
                if (atomic_group) {
                    ++stats.atomic_attempts;
                    bool accepted = false;
                    if (config.atomic_capacity_candidates > 1) {
                        std::sort(atomic_candidates.begin(), atomic_candidates.end(),
                                  [](const AtomicCandidate& a, const AtomicCandidate& b) {
                            if (a.individual_score != b.individual_score) {
                                return a.individual_score < b.individual_score;
                            }
                            return a.bin < b.bin;
                        });
                        if (static_cast<int>(atomic_candidates.size()) >
                            config.atomic_capacity_candidates) {
                            atomic_candidates.resize(config.atomic_capacity_candidates);
                        }
                        bool found_bid = false;
                        bool first_bid_failed = false;
                        Real best_bid_hpwl = std::numeric_limits<Real>::infinity();
                        AtomicCandidate best_bid;
                        for (std::size_t bid_index = 0;
                             bid_index < atomic_candidates.size(); ++bid_index) {
                            const AtomicCandidate& candidate = atomic_candidates[bid_index];
                            const Real dx = candidate.x - node.x;
                            const Real dy = candidate.y - node.y;
                            const Real maximum_group_delta = config.atomic_min_gain_ratio *
                                candidate.first_overflow_delta;
                            const AtomicMoveResult trial = apply_atomic_group_move(
                                db, density, *atomic_members, dx, dy,
                                maximum_group_delta, tolerance, true, false, occupancy);
                            ++stats.atomic_bid_trials;
                            if (bid_index == 0 && trial != AtomicMoveResult::Accepted) {
                                first_bid_failed = true;
                            }
                            if (trial == AtomicMoveResult::CapacityRejected) {
                                ++stats.atomic_capacity_rejections;
                            }
                            if (trial != AtomicMoveResult::Accepted) continue;
                            ++stats.atomic_feasible_bids;
                            const Real hpwl_delta = group_incident_hpwl_delta(
                                db, *atomic_members, dx, dy, node_nets);
                            if (!found_bid || hpwl_delta < best_bid_hpwl) {
                                found_bid = true;
                                best_bid_hpwl = hpwl_delta;
                                best_bid = candidate;
                            }
                        }
                        if (found_bid) {
                            best_x = best_bid.x;
                            best_y = best_bid.y;
                            best_bin = best_bid.bin;
                            const Real dx = best_x - node.x;
                            const Real dy = best_y - node.y;
                            const Real maximum_group_delta = config.atomic_min_gain_ratio *
                                best_bid.first_overflow_delta;
                            const AtomicMoveResult committed = apply_atomic_group_move(
                                db, density, *atomic_members, dx, dy,
                                maximum_group_delta, tolerance, true, true, occupancy);
                            if (committed != AtomicMoveResult::Accepted) {
                                throw std::logic_error("atomic capacity bid changed before commit");
                            }
                            accepted = true;
                            if (first_bid_failed) ++stats.atomic_rescued_groups;
                        }
                    } else {
                        const Real dx = best_x - node.x;
                        const Real dy = best_y - node.y;
                        const Real maximum_group_delta =
                            config.atomic_min_gain_ratio * best_delta.overflow_area;
                        const AtomicMoveResult atomic_result = apply_atomic_group_move(
                            db, density, *atomic_members, dx, dy,
                            maximum_group_delta, tolerance,
                            config.atomic_componentwise_capacity, true, occupancy);
                        accepted = atomic_result == AtomicMoveResult::Accepted;
                        if (atomic_result == AtomicMoveResult::CapacityRejected) {
                            ++stats.atomic_capacity_rejections;
                        }
                    }
                    if (!accepted && config.atomic_fallback_individual) {
                        atomic_group_fallback[group_id] = 1;
                        ++stats.atomic_fallbacks;
                        apply_move(node, best_x, best_y, best_delta, occupancy);
                        moved[id] = 1;
                        ++stats.moves;
                        ++round_moves;
                        if (config.group_destination_radius > 0 &&
                            group_destination[group_id] < 0) {
                            group_destination[group_id] = best_bin;
                        }
                    }
                    for (const Destination& destination : reusable) {
                        if (occupancy[destination.bin] < capacity - tolerance) {
                            heap.push(destination);
                        }
                    }
                    if (!accepted) continue;
                    if (config.group_destination_radius > 0 &&
                        group_destination[group_id] < 0) {
                        group_destination[group_id] = best_bin;
                    }
                    for (int member : *atomic_members) moved[member] = 1;
                    const int accepted_nodes = static_cast<int>(atomic_members->size());
                    stats.moves += accepted_nodes;
                    round_moves += accepted_nodes;
                    ++stats.atomic_groups;
                    continue;
                }
                if (rigid_groups && group_members[group_id] > 1 &&
                    !group_has_displacement[group_id]) {
                    group_dx[group_id] = best_x - node.x;
                    group_dy[group_id] = best_y - node.y;
                    group_has_displacement[group_id] = 1;
                }
                apply_move(node, best_x, best_y, best_delta, occupancy);
                for (const Destination& destination : reusable) {
                    if (occupancy[destination.bin] < capacity - tolerance) heap.push(destination);
                }
                if (config.group_destination_radius > 0 &&
                    group_members[group_id] > 1 && group_destination[group_id] < 0) {
                    group_destination[group_id] = best_bin;
                }
                moved[id] = 1;
                ++stats.moves;
                ++round_moves;
            }
        }
        if (config.global_hpwl_budget >= 0.0) {
            const Real round_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
            if (round_hpwl > round_start_hpwl * (1.0 + config.global_hpwl_budget)) {
                for (const auto& item : round_positions) {
                    db.nodes[item.first].x = item.second.first;
                    db.nodes[item.first].y = item.second.second;
                }
                occupancy = round_occupancy;
                stats.moves = moves_before_round;
                round_moves = 0;
            }
        }
        ++stats.rounds;
        if (round_moves == 0) break;
    }

    const DensityMetrics final_density = density.evaluate(0.0, 1.0, nullptr, nullptr);
    stats.final_overflow = final_density.overflow;
    stats.final_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    stats.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return stats;
}

}  // namespace ea
