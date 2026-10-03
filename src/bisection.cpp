#include "epsilon_active/bisection.hpp"

#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ea {
namespace {

struct Region {
    int x0 = 0;
    int x1 = 0;
    int y0 = 0;
    int y1 = 0;
    int depth = 0;
    std::vector<int> nodes;
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

Real bin_capacity(const ExactOverlapDensity& density, int bin) {
    return std::max<Real>(density.target_density() * density.bin_area() -
                          density.fixed_occupancy()[bin], 0.0);
}

Real region_capacity(const ExactOverlapDensity& density, const Region& region) {
    Real result = 0.0;
    for (int by = region.y0; by < region.y1; ++by) {
        for (int bx = region.x0; bx < region.x1; ++bx) {
            result += bin_capacity(density, by * density.bins_x() + bx);
        }
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

class Bisector {
public:
    Bisector(Database& db, ExactOverlapDensity& density,
             const BisectionConfig& config, BisectionStats& stats)
        : db_(db), density_(density), config_(config), stats_(stats),
          node_nets_(build_node_nets(db)), membership_(db.nodes.size(), 0),
          side_(db.nodes.size(), -1), score_(db.nodes.size(), 0.0),
          net_left_(db.nets.size(), 0), net_right_(db.nets.size(), 0),
          leaf_of_(db.nodes.size(), -1),
          remaining_capacity_(static_cast<std::size_t>(density.bins_x()) *
                              density.bins_y(), 0.0),
          atomic_group_of_(db.nodes.size(), -1),
          unit_of_node_(db.nodes.size(), -1) {
        for (int bin = 0; bin < static_cast<int>(remaining_capacity_.size());
             ++bin) {
            remaining_capacity_[bin] = bin_capacity(density_, bin);
        }
        build_atomic_groups();
    }

    void run() {
        if (config_.axis_rank_map) apply_axis_rank_map();
        if (config_.curve_rank_map) apply_curve_rank_map();
        if (config_.axis_map_only || config_.curve_map_only) return;
        Region root{0, density_.bins_x(), 0, density_.bins_y(), 0,
                    db_.movable_ids};
        std::vector<Region> stack;
        stack.push_back(std::move(root));
        std::vector<Region> leaves;
        while (!stack.empty()) {
            Region region = std::move(stack.back());
            stack.pop_back();
            stats_.max_depth = std::max(stats_.max_depth, region.depth);
            const int width = region.x1 - region.x0;
            const int height = region.y1 - region.y0;
            if (region.nodes.size() < 2 ||
                (config_.max_depth >= 0 && region.depth >= config_.max_depth) ||
                (width <= config_.leaf_bins && height <= config_.leaf_bins)) {
                leaves.push_back(std::move(region));
                continue;
            }
            auto [first, second] = split(std::move(region));
            if (first.nodes.empty()) {
                stack.push_back(std::move(second));
                continue;
            }
            if (second.nodes.empty()) {
                stack.push_back(std::move(first));
                continue;
            }
            ++stats_.splits;
            stack.push_back(std::move(second));
            stack.push_back(std::move(first));
        }
        stats_.leaves = static_cast<int>(leaves.size());
        for (int leaf = 0; leaf < static_cast<int>(leaves.size()); ++leaf) {
            for (int id : leaves[leaf].nodes) leaf_of_[id] = leaf;
            place_leaf(leaves[leaf]);
            const Real capacity = region_capacity(density_, leaves[leaf]);
            Real area = 0.0;
            for (int id : leaves[leaf].nodes) area += db_.nodes[id].area();
            stats_.maximum_capacity_error = std::max(
                stats_.maximum_capacity_error,
                std::max<Real>(area - capacity, 0.0) /
                    std::max<Real>(capacity, 1.0));
        }
        measure_cut_weight();
    }

private:
    Database& db_;
    ExactOverlapDensity& density_;
    const BisectionConfig& config_;
    BisectionStats& stats_;
    std::vector<std::vector<int>> node_nets_;
    std::vector<int> membership_;
    std::vector<signed char> side_;
    std::vector<Real> score_;
    std::vector<int> net_left_;
    std::vector<int> net_right_;
    std::vector<int> leaf_of_;
    std::vector<Real> remaining_capacity_;
    std::vector<int> atomic_group_of_;
    std::vector<std::vector<int>> atomic_groups_;
    std::vector<int> unit_of_node_;
    int stamp_ = 0;

    void build_atomic_groups() {
        if (!config_.atomic_coarsening) return;
        std::vector<int> parent(db_.nodes.size(), -1);
        std::vector<int> group_size(db_.nodes.size(), 0);
        for (int id : db_.movable_ids) {
            parent[id] = id;
            group_size[id] = 1;
        }
        const auto find_root = [&](int id, const auto& self) -> int {
            if (parent[id] == id) return id;
            parent[id] = self(parent[id], self);
            return parent[id];
        };
        for (int round = 0; round < config_.atomic_rounds; ++round) {
            std::vector<int> best(db_.nodes.size(), -1);
            std::vector<Real> best_score(
                db_.nodes.size(), -std::numeric_limits<Real>::infinity());
            // Accumulate all shared-net contributions for a root pair before
            // selecting mutual heavy-edge matches.  The previous implementation
            // kept only the largest single net contribution, which made the
            // coarsening almost ineffective on multi-net connected cells.
            struct PairContribution {
                int a = -1;
                int b = -1;
                Real weight = 0.0;
            };
            std::vector<PairContribution> contributions;
            contributions.reserve(db_.nets.size() * 2);
            for (const Net& net : db_.nets) {
                if (net.pin_count < 2 ||
                    net.pin_count >
                        static_cast<std::size_t>(config_.atomic_degree_limit)) {
                    continue;
                }
                std::vector<int> roots;
                roots.reserve(net.pin_count);
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const int id = db_.pins[p].node;
                    if (db_.nodes[id].fixed) continue;
                    roots.push_back(find_root(id, find_root));
                }
                std::sort(roots.begin(), roots.end());
                roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
                const Real edge_score = net.weight /
                    std::max<Real>(static_cast<Real>(roots.size()) - 1.0, 1.0);
                for (std::size_t i = 0; i < roots.size(); ++i) {
                    for (std::size_t j = i + 1; j < roots.size(); ++j) {
                        const int a = roots[i];
                        const int b = roots[j];
                        if (group_size[a] + group_size[b] >
                            config_.atomic_max_nodes) continue;
                        contributions.push_back({a, b, edge_score});
                    }
                }
            }
            std::sort(contributions.begin(), contributions.end(),
                      [](const PairContribution& lhs,
                         const PairContribution& rhs) {
                if (lhs.a != rhs.a) return lhs.a < rhs.a;
                return lhs.b < rhs.b;
            });
            for (std::size_t i = 0; i < contributions.size();) {
                const int a = contributions[i].a;
                const int b = contributions[i].b;
                Real score = 0.0;
                std::size_t j = i;
                while (j < contributions.size() &&
                       contributions[j].a == a && contributions[j].b == b) {
                    score += contributions[j].weight;
                    ++j;
                }
                if (score > best_score[a] ||
                    (score == best_score[a] && b < best[a])) {
                    best_score[a] = score;
                    best[a] = b;
                }
                if (score > best_score[b] ||
                    (score == best_score[b] && a < best[b])) {
                    best_score[b] = score;
                    best[b] = a;
                }
                i = j;
            }
            int merges = 0;
            for (int id : db_.movable_ids) {
                const int root = find_root(id, find_root);
                if (root != id || best[root] < 0) continue;
                const int other = find_root(best[root], find_root);
                if (other == root || best[other] != root ||
                    group_size[root] + group_size[other] >
                        config_.atomic_max_nodes) continue;
                const int keep = std::min(root, other);
                const int remove = std::max(root, other);
                parent[remove] = keep;
                group_size[keep] += group_size[remove];
                group_size[remove] = 0;
                ++merges;
                ++stats_.atomic_merges;
            }
            if (merges == 0) break;
        }
        std::vector<int> root_to_group(db_.nodes.size(), -1);
        for (int id : db_.movable_ids) {
            const int root = find_root(id, find_root);
            if (root_to_group[root] < 0) {
                root_to_group[root] = static_cast<int>(atomic_groups_.size());
                atomic_groups_.push_back({});
            }
            const int group = root_to_group[root];
            atomic_group_of_[id] = group;
            atomic_groups_[group].push_back(id);
        }
        stats_.atomic_groups = static_cast<int>(atomic_groups_.size());
    }

    void apply_axis_rank_map() {
        const auto map_axis = [&](bool vertical) {
            const int count = vertical ? density_.bins_x() : density_.bins_y();
            const Real pitch = vertical ? density_.bin_width() : density_.bin_height();
            std::vector<Real> cumulative(count, 0.0);
            Real total_capacity = 0.0;
            for (int index = 0; index < count; ++index) {
                Real capacity = 0.0;
                if (vertical) {
                    for (int by = 0; by < density_.bins_y(); ++by) {
                        capacity += bin_capacity(density_,
                            by * density_.bins_x() + index);
                    }
                } else {
                    for (int bx = 0; bx < density_.bins_x(); ++bx) {
                        capacity += bin_capacity(density_,
                            index * density_.bins_x() + bx);
                    }
                }
                total_capacity += std::max<Real>(capacity, 1.0e-12);
                cumulative[index] = total_capacity;
            }
            if (total_capacity <= 0.0 || db_.movable_ids.empty()) return;
            std::vector<int> order = db_.movable_ids;
            std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
                const Real left = vertical ? db_.nodes[lhs].x : db_.nodes[lhs].y;
                const Real right = vertical ? db_.nodes[rhs].x : db_.nodes[rhs].y;
                if (left != right) return left < right;
                return lhs < rhs;
            });
            Real total_area = 0.0;
            for (int id : order) total_area += db_.nodes[id].area();
            Real prefix = 0.0;
            for (int id : order) {
                const Real quantile = total_area > 0.0
                    ? (prefix + 0.5 * db_.nodes[id].area()) / total_area * total_capacity
                    : 0.0;
                const auto found = std::lower_bound(cumulative.begin(), cumulative.end(),
                                                    quantile);
                const int index = std::min<int>(
                    static_cast<int>(std::distance(cumulative.begin(), found)), count - 1);
                const Real coordinate = (vertical ? db_.xl : db_.yl) +
                    (index + 0.5) * pitch;
                Node& node = db_.nodes[id];
                if (vertical) {
                    node.x = std::clamp(coordinate, db_.xl + 0.5 * node.width,
                                        db_.xh - 0.5 * node.width);
                } else {
                    node.y = std::clamp(coordinate, db_.yl + 0.5 * node.height,
                                        db_.yh - 0.5 * node.height);
                }
                prefix += node.area();
            }
        };
        map_axis(true);
        map_axis(false);
    }

    void apply_curve_rank_map() {
        const int bins_x = density_.bins_x();
        const int bins_y = density_.bins_y();
        const bool hilbert_grid = bins_x == bins_y &&
            (bins_x & (bins_x - 1)) == 0;
        std::vector<int> bins;
        bins.reserve(bins_x * bins_y);
        for (int by = 0; by < bins_y; ++by) {
            for (int bx = 0; bx < bins_x; ++bx) bins.push_back(by * bins_x + bx);
        }
        std::sort(bins.begin(), bins.end(), [&](int lhs, int rhs) {
            if (!hilbert_grid) return lhs < rhs;
            const auto left = hilbert_index(bins_x, lhs % bins_x, lhs / bins_x);
            const auto right = hilbert_index(bins_x, rhs % bins_x, rhs / bins_x);
            if (left != right) return left < right;
            return lhs < rhs;
        });
        std::vector<Real> cumulative(bins.size(), 0.0);
        Real total_capacity = 0.0;
        for (std::size_t i = 0; i < bins.size(); ++i) {
            total_capacity += std::max<Real>(bin_capacity(density_, bins[i]), 1.0e-12);
            cumulative[i] = total_capacity;
        }
        if (db_.movable_ids.empty()) return;
        std::vector<int> order = db_.movable_ids;
        const auto curve_key = [&](int id) {
            const int bx = std::clamp(static_cast<int>(
                (db_.nodes[id].x - db_.xl) / density_.bin_width()), 0, bins_x - 1);
            const int by = std::clamp(static_cast<int>(
                (db_.nodes[id].y - db_.yl) / density_.bin_height()), 0, bins_y - 1);
            if (hilbert_grid) return hilbert_index(bins_x, bx, by);
            return static_cast<std::uint64_t>(by * bins_x + bx);
        };
        std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
            const auto left = curve_key(lhs);
            const auto right = curve_key(rhs);
            if (left != right) return left < right;
            return lhs < rhs;
        });
        Real total_area = 0.0;
        for (int id : order) total_area += db_.nodes[id].area();
        Real prefix = 0.0;
        for (int id : order) {
            const Real quantile = total_area > 0.0
                ? (prefix + 0.5 * db_.nodes[id].area()) / total_area * total_capacity
                : 0.0;
            const auto found = std::lower_bound(cumulative.begin(), cumulative.end(),
                                                quantile);
            const int position = std::min<int>(
                static_cast<int>(std::distance(cumulative.begin(), found)),
                static_cast<int>(bins.size()) - 1);
            const int bin = bins[position];
            const int bx = bin % bins_x;
            const int by = bin / bins_x;
            Node& node = db_.nodes[id];
            node.x = std::clamp(db_.xl + (bx + 0.5) * density_.bin_width(),
                                db_.xl + 0.5 * node.width,
                                db_.xh - 0.5 * node.width);
            node.y = std::clamp(db_.yl + (by + 0.5) * density_.bin_height(),
                                db_.yl + 0.5 * node.height,
                                db_.yh - 0.5 * node.height);
            prefix += node.area();
        }
    }

    std::pair<Region, Region> split(Region region) {
        const bool vertical = (region.x1 - region.x0) >=
                              (region.y1 - region.y0) &&
                              region.x1 - region.x0 > 1;
        const int low = vertical ? region.x0 : region.y0;
        const int high = vertical ? region.x1 : region.y1;
        const Real total_capacity = region_capacity(density_, region);
        Real prefix = 0.0;
        Real best_error = std::numeric_limits<Real>::infinity();
        int coordinate = low + 1;
        for (int cut = low + 1; cut < high; ++cut) {
            const int slice = cut - 1;
            if (vertical) {
                for (int by = region.y0; by < region.y1; ++by) {
                    prefix += bin_capacity(
                        density_, by * density_.bins_x() + slice);
                }
            } else {
                for (int bx = region.x0; bx < region.x1; ++bx) {
                    prefix += bin_capacity(
                        density_, slice * density_.bins_x() + bx);
                }
            }
            const Real error = std::abs(prefix - 0.5 * total_capacity);
            if (error < best_error) {
                best_error = error;
                coordinate = cut;
            }
        }
        Region first = region;
        Region second = region;
        first.depth = second.depth = region.depth + 1;
        first.nodes.clear();
        second.nodes.clear();
        if (vertical) {
            first.x1 = coordinate;
            second.x0 = coordinate;
        } else {
            first.y1 = coordinate;
            second.y0 = coordinate;
        }
        const Real first_capacity = region_capacity(density_, first);
        const Real second_capacity = region_capacity(density_, second);
        partition_nodes(region.nodes, vertical, coordinate, region.depth,
                        first_capacity, second_capacity,
                        first.nodes, second.nodes);
        return {std::move(first), std::move(second)};
    }

    bool partition_atomic(const std::vector<int>& nodes, bool vertical,
                          int coordinate, Real target_first, Real total_area,
                          Real first_capacity, Real second_capacity,
                          Real maximum_node_area, std::vector<int>& first,
                          std::vector<int>& second) {
        struct Unit {
            std::vector<int> nodes;
            Real area = 0.0;
            Real coordinate = 0.0;
            Real seed_score = 0.0;
            signed char side = -1;
        };
        std::vector<Unit> grouped;
        std::vector<int> group_to_unit(atomic_groups_.size(), -1);
        for (int id : nodes) {
            const int group = atomic_group_of_[id];
            int unit = group >= 0 ? group_to_unit[group] : -1;
            if (unit < 0) {
                unit = static_cast<int>(grouped.size());
                grouped.push_back({});
                if (group >= 0) group_to_unit[group] = unit;
            }
            grouped[unit].nodes.push_back(id);
        }

        std::vector<Unit> units;
        for (Unit& group : grouped) {
            Real area = 0.0;
            for (int id : group.nodes) area += db_.nodes[id].area();
            if (area > std::max(first_capacity, second_capacity) + 1.0e-9) {
                for (int id : group.nodes) {
                    Unit singleton;
                    singleton.nodes.push_back(id);
                    units.push_back(std::move(singleton));
                }
            } else {
                units.push_back(std::move(group));
            }
        }
        if (units.size() < 2) {
            return false;
        }
        if (target_first <= 1.0e-12 ||
            total_area - target_first <= 1.0e-12) {
            return false;
        }

        Real maximum_unit_area = maximum_node_area;
        for (int unit = 0; unit < static_cast<int>(units.size()); ++unit) {
            Real weighted_coordinate = 0.0;
            for (int id : units[unit].nodes) {
                const Real area = db_.nodes[id].area();
                units[unit].area += area;
                weighted_coordinate += area * (vertical
                    ? (db_.nodes[id].x - db_.xl) / density_.bin_width()
                    : (db_.nodes[id].y - db_.yl) / density_.bin_height());
                units[unit].seed_score += score_[id];
                unit_of_node_[id] = unit;
            }
            units[unit].coordinate = weighted_coordinate /
                std::max<Real>(units[unit].area, 1.0e-30);
            maximum_unit_area = std::max(maximum_unit_area, units[unit].area);
        }

        std::vector<int> order(units.size());
        std::iota(order.begin(), order.end(), 0);
        if (config_.position_seeded) {
            std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
                if (units[lhs].coordinate != units[rhs].coordinate) {
                    return units[lhs].coordinate < units[rhs].coordinate;
                }
                return lhs < rhs;
            });
        } else {
            std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
                if (units[lhs].seed_score != units[rhs].seed_score) {
                    return units[lhs].seed_score > units[rhs].seed_score;
                }
                return lhs < rhs;
            });
        }

        Real first_area = 0.0;
        for (std::size_t position = 0; position < order.size(); ++position) {
            const int unit = order[position];
            const bool leave_one_for_second = position + 1 == order.size();
            const bool need_first = first_area < target_first &&
                !leave_one_for_second;
            if (need_first &&
                first_area + units[unit].area <= first_capacity + 1.0e-9) {
                units[unit].side = 0;
                first_area += units[unit].area;
            } else {
                units[unit].side = 1;
            }
        }
        if (first_area <= 0.0) {
            const int unit = order.front();
            units[unit].side = 0;
            first_area += units[unit].area;
        }

        auto move_unit = [&](int unit, signed char destination) {
            if (units[unit].side == destination) return;
            if (destination == 0) first_area += units[unit].area;
            else first_area -= units[unit].area;
            units[unit].side = destination;
        };
        while (total_area - first_area > second_capacity + 1.0e-9) {
            int best = -1;
            for (int unit : order) {
                if (units[unit].side != 1 ||
                    first_area + units[unit].area > first_capacity + 1.0e-9)
                    continue;
                if (best < 0 || units[unit].area < units[best].area) best = unit;
            }
            if (best < 0) break;
            move_unit(best, 0);
        }
        while (first_area > first_capacity + 1.0e-9) {
            int best = -1;
            for (int unit : order) {
                if (units[unit].side != 0 ||
                    total_area - first_area + units[unit].area >
                        second_capacity + 1.0e-9) continue;
                if (best < 0 || units[unit].area < units[best].area) best = unit;
            }
            if (best < 0) break;
            move_unit(best, 1);
        }

        std::vector<int> touched;
        for (int id : nodes) {
            side_[id] = units[unit_of_node_[id]].side;
            for (int net_id : node_nets_[id]) {
                if (net_left_[net_id] == 0 && net_right_[net_id] == 0) {
                    touched.push_back(net_id);
                }
                if (side_[id] == 0) ++net_left_[net_id];
                else ++net_right_[net_id];
            }
        }
        for (int net_id : touched) {
            const Net& net = db_.nets[net_id];
            for (std::size_t p = net.pin_begin;
                 p < net.pin_begin + net.pin_count; ++p) {
                const Node& node = db_.nodes[db_.pins[p].node];
                if (!node.fixed) continue;
                const Real value = vertical
                    ? (node.x - db_.xl) / density_.bin_width()
                    : (node.y - db_.yl) / density_.bin_height();
                if (value < coordinate) ++net_left_[net_id];
                else ++net_right_[net_id];
            }
        }
        const Real tolerance = std::max(
            config_.balance_tolerance * total_area, maximum_unit_area);
        const auto unit_gain = [&](int unit) {
            std::vector<int> nets;
            for (int id : units[unit].nodes) {
                nets.insert(nets.end(), node_nets_[id].begin(),
                            node_nets_[id].end());
            }
            std::sort(nets.begin(), nets.end());
            nets.erase(std::unique(nets.begin(), nets.end()), nets.end());
            Real gain = 0.0;
            for (int net_id : nets) {
                int unit_pins = 0;
                const Net& net = db_.nets[net_id];
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const int id = db_.pins[p].node;
                    if (!db_.nodes[id].fixed && membership_[id] == stamp_ &&
                        unit_of_node_[id] == unit) ++unit_pins;
                }
                int left = net_left_[net_id];
                int right = net_right_[net_id];
                const bool old_cut = left > 0 && right > 0;
                if (units[unit].side == 0) {
                    left -= unit_pins;
                    right += unit_pins;
                } else {
                    right -= unit_pins;
                    left += unit_pins;
                }
                const bool new_cut = left > 0 && right > 0;
                gain += net.weight *
                    (static_cast<Real>(old_cut) - static_cast<Real>(new_cut));
            }
            return gain;
        };

        for (int pass = 0; pass < config_.fm_passes; ++pass) {
            std::vector<std::pair<Real, int>> gains;
            gains.reserve(units.size());
            for (int unit = 0; unit < static_cast<int>(units.size()); ++unit) {
                gains.push_back({unit_gain(unit), unit});
            }
            std::sort(gains.begin(), gains.end(), [](const auto& lhs,
                                                     const auto& rhs) {
                if (lhs.first != rhs.first) return lhs.first > rhs.first;
                return lhs.second < rhs.second;
            });
            int moved = 0;
            for (const auto& [stale_gain, unit] : gains) {
                (void)stale_gain;
                const Real gain = unit_gain(unit);
                if (gain <= 1.0e-12) continue;
                const Real new_first = first_area +
                    (units[unit].side == 0
                        ? -units[unit].area : units[unit].area);
                const Real new_second = total_area - new_first;
                if (new_first > first_capacity + 1.0e-9 ||
                    new_second > second_capacity + 1.0e-9 ||
                    new_first < target_first - tolerance ||
                    new_first > target_first + tolerance) continue;
                for (int id : units[unit].nodes) {
                    for (int net_id : node_nets_[id]) {
                        if (units[unit].side == 0) {
                            --net_left_[net_id];
                            ++net_right_[net_id];
                        } else {
                            --net_right_[net_id];
                            ++net_left_[net_id];
                        }
                    }
                }
                move_unit(unit, units[unit].side == 0 ? 1 : 0);
                for (int id : units[unit].nodes) {
                    side_[id] = units[unit].side;
                }
                ++moved;
                ++stats_.atomic_unit_moves;
                stats_.fm_moves += static_cast<int>(units[unit].nodes.size());
            }
            if (moved == 0) break;
        }
        for (int net_id : touched) {
            net_left_[net_id] = 0;
            net_right_[net_id] = 0;
        }
        for (int id : nodes) {
            (side_[id] == 0 ? first : second).push_back(id);
        }
        return !first.empty() && !second.empty();
    }

    void partition_nodes(const std::vector<int>& nodes, bool vertical,
                         int coordinate, int depth, Real first_capacity,
                         Real second_capacity, std::vector<int>& first,
                         std::vector<int>& second) {
        ++stamp_;
        Real total_area = 0.0;
        Real maximum_area = 0.0;
        for (int id : nodes) {
            membership_[id] = stamp_;
            side_[id] = -1;
            score_[id] = 0.0;
            total_area += db_.nodes[id].area();
            maximum_area = std::max(maximum_area, db_.nodes[id].area());
        }
        const Real capacity_sum = first_capacity + second_capacity;
        const Real target_first = capacity_sum > 0.0
            ? total_area * first_capacity / capacity_sum
            : 0.5 * total_area;
        for (int id : nodes) {
            for (int net_id : node_nets_[id]) {
                const Net& net = db_.nets[net_id];
                bool fixed_first = false;
                bool fixed_second = false;
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const Node& node = db_.nodes[db_.pins[p].node];
                    if (!node.fixed) continue;
                    const Real value = vertical
                        ? (node.x - db_.xl) / density_.bin_width()
                        : (node.y - db_.yl) / density_.bin_height();
                    fixed_first = fixed_first || value < coordinate;
                    fixed_second = fixed_second || value >= coordinate;
                }
                if (fixed_first != fixed_second) {
                    score_[id] += fixed_first ? net.weight : -net.weight;
                }
            }
        }

        if (config_.atomic_coarsening &&
            depth < config_.atomic_release_depth) {
            if (partition_atomic(nodes, vertical, coordinate, target_first,
                                 total_area, first_capacity, second_capacity,
                                 maximum_area, first, second)) {
                return;
            }
        }

        if (config_.position_seeded) {
            std::vector<int> first_candidates;
            std::vector<int> second_candidates;
            first_candidates.reserve(nodes.size());
            second_candidates.reserve(nodes.size());
            for (int id : nodes) {
                const Real value = vertical
                    ? (db_.nodes[id].x - db_.xl) / density_.bin_width()
                    : (db_.nodes[id].y - db_.yl) / density_.bin_height();
                if (value < coordinate) {
                    side_[id] = 0;
                    first_candidates.push_back(id);
                } else {
                    side_[id] = 1;
                    second_candidates.push_back(id);
                }
            }
            auto by_boundary = [&](int lhs, int rhs) {
                const Real left = vertical
                    ? (db_.nodes[lhs].x - db_.xl) / density_.bin_width()
                    : (db_.nodes[lhs].y - db_.yl) / density_.bin_height();
                const Real right = vertical
                    ? (db_.nodes[rhs].x - db_.xl) / density_.bin_width()
                    : (db_.nodes[rhs].y - db_.yl) / density_.bin_height();
                const Real left_distance = std::abs(left - coordinate);
                const Real right_distance = std::abs(right - coordinate);
                if (left_distance != right_distance) return left_distance < right_distance;
                return lhs < rhs;
            };
            std::sort(first_candidates.begin(), first_candidates.end(), by_boundary);
            std::sort(second_candidates.begin(), second_candidates.end(), by_boundary);
            Real first_area = 0.0;
            for (int id : nodes) if (side_[id] == 0) first_area += db_.nodes[id].area();
            if (config_.surplus_only) {
                // Keep the position-induced side assignment.  Repair only the
                // side whose exact capacity is exceeded; unlike target
                // balancing this does not move cells merely to equalize a cut.
                auto move_surplus_first = [&]() {
                    for (int id : first_candidates) {
                        if (first_area <= first_capacity + 1.0e-9) break;
                        side_[id] = 1;
                        first_area -= db_.nodes[id].area();
                    }
                };
                auto move_surplus_second = [&]() {
                    for (int id : second_candidates) {
                        if (total_area - first_area <= second_capacity + 1.0e-9) break;
                        side_[id] = 0;
                        first_area += db_.nodes[id].area();
                    }
                };
                if (first_area > first_capacity + 1.0e-9) {
                    move_surplus_first();
                }
                if (total_area - first_area > second_capacity + 1.0e-9) {
                    move_surplus_second();
                }
                for (int id : nodes) {
                    (side_[id] == 0 ? first : second).push_back(id);
                }
                return;
            }
            const auto move_first_to_second = [&]() {
                for (int id : first_candidates) {
                    if (first_candidates.size() <= 1 || first_area <= target_first) break;
                    side_[id] = 1;
                    first_area -= db_.nodes[id].area();
                    if (first_area <= target_first) break;
                }
            };
            const auto move_second_to_first = [&]() {
                for (int id : second_candidates) {
                    if (second_candidates.size() <= 1 || first_area >= target_first) break;
                    side_[id] = 0;
                    first_area += db_.nodes[id].area();
                    if (first_area >= target_first) break;
                }
            };
            if (first_area > target_first) move_first_to_second();
            else move_second_to_first();
            if (std::all_of(nodes.begin(), nodes.end(),
                            [&](int id) { return side_[id] == 0; })) {
                side_[nodes.back()] = 1;
            } else if (std::all_of(nodes.begin(), nodes.end(),
                                   [&](int id) { return side_[id] == 1; })) {
                side_[nodes.front()] = 0;
            }
            for (int id : nodes) {
                (side_[id] == 0 ? first : second).push_back(id);
            }
            refine_fm(nodes, vertical, coordinate, target_first, total_area,
                      first_capacity, second_capacity,
                      std::max(config_.balance_tolerance * total_area,
                               maximum_area));
            first.clear();
            second.clear();
            for (int id : nodes) {
                (side_[id] == 0 ? first : second).push_back(id);
            }
            return;
        }

        using Candidate = std::pair<Real, int>;
        std::priority_queue<Candidate> queue;
        for (int id : nodes) queue.push({score_[id], -id});
        Real first_area = 0.0;
        while (first.size() + 1 < nodes.size() && first_area < target_first) {
            int id = -1;
            while (!queue.empty()) {
                const auto [queued_score, negated] = queue.top();
                queue.pop();
                const int candidate = -negated;
                if (side_[candidate] < 0 &&
                    std::abs(queued_score - score_[candidate]) <= 1.0e-12) {
                    id = candidate;
                    break;
                }
            }
            if (id < 0) {
                const auto found = std::find_if(nodes.begin(), nodes.end(),
                    [&](int candidate) { return side_[candidate] < 0; });
                if (found == nodes.end()) break;
                id = *found;
            }
            side_[id] = 0;
            first.push_back(id);
            first_area += db_.nodes[id].area();
            for (int net_id : node_nets_[id]) {
                const Net& net = db_.nets[net_id];
                if (net.pin_count > static_cast<std::size_t>(config_.degree_limit)) {
                    continue;
                }
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const int neighbor = db_.pins[p].node;
                    if (membership_[neighbor] == stamp_ && side_[neighbor] < 0) {
                        score_[neighbor] += net.weight;
                        queue.push({score_[neighbor], -neighbor});
                    }
                }
            }
        }
        for (int id : nodes) {
            if (side_[id] < 0) {
                side_[id] = 1;
                second.push_back(id);
            }
        }
        refine_fm(nodes, vertical, coordinate, target_first, total_area,
                  first_capacity, second_capacity,
                  std::max(config_.balance_tolerance * total_area,
                           maximum_area));
        first.clear();
        second.clear();
        for (int id : nodes) {
            (side_[id] == 0 ? first : second).push_back(id);
        }
    }

    void refine_fm(const std::vector<int>& nodes, bool vertical, int coordinate,
                   Real target_first, Real total_area, Real first_capacity,
                   Real second_capacity, Real tolerance) {
        std::vector<int> touched;
        for (int id : nodes) {
            for (int net_id : node_nets_[id]) {
                if (net_left_[net_id] == 0 && net_right_[net_id] == 0) {
                    touched.push_back(net_id);
                }
                if (side_[id] == 0) ++net_left_[net_id];
                else ++net_right_[net_id];
            }
        }
        for (int net_id : touched) {
            const Net& net = db_.nets[net_id];
            for (std::size_t p = net.pin_begin;
                 p < net.pin_begin + net.pin_count; ++p) {
                const Node& node = db_.nodes[db_.pins[p].node];
                if (!node.fixed) continue;
                const Real value = vertical
                    ? (node.x - db_.xl) / density_.bin_width()
                    : (node.y - db_.yl) / density_.bin_height();
                if (value < coordinate) ++net_left_[net_id];
                else ++net_right_[net_id];
            }
        }
        Real first_area = 0.0;
        for (int id : nodes) if (side_[id] == 0) first_area += db_.nodes[id].area();
        const Real target_second = total_area - target_first;
        for (int pass = 0; pass < config_.fm_passes; ++pass) {
            std::vector<std::pair<Real, int>> gains;
            gains.reserve(nodes.size());
            for (int id : nodes) gains.push_back({fm_gain(id), id});
            std::sort(gains.begin(), gains.end(), [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second < b.second;
            });
            int moved = 0;
            for (const auto& [stale_gain, id] : gains) {
                (void)stale_gain;
                const Real gain = fm_gain(id);
                if (gain <= 1.0e-12) continue;
                const Real area = db_.nodes[id].area();
                const Real new_first = first_area + (side_[id] == 0 ? -area : area);
                const Real new_second = total_area - new_first;
                if (new_first > first_capacity + 1.0e-9 ||
                    new_second > second_capacity + 1.0e-9 ||
                    new_first < target_first - tolerance ||
                    new_first > target_first + tolerance ||
                    new_second < target_second - tolerance ||
                    new_second > target_second + tolerance) {
                    continue;
                }
                for (int net_id : node_nets_[id]) {
                    if (side_[id] == 0) {
                        --net_left_[net_id];
                        ++net_right_[net_id];
                    } else {
                        --net_right_[net_id];
                        ++net_left_[net_id];
                    }
                }
                side_[id] = side_[id] == 0 ? 1 : 0;
                first_area = new_first;
                ++moved;
                ++stats_.fm_moves;
            }
            if (moved == 0) break;
        }
        for (int net_id : touched) {
            net_left_[net_id] = 0;
            net_right_[net_id] = 0;
        }
    }

    Real fm_gain(int id) const {
        Real gain = 0.0;
        for (int net_id : node_nets_[id]) {
            const bool old_cut = net_left_[net_id] > 0 && net_right_[net_id] > 0;
            int left = net_left_[net_id];
            int right = net_right_[net_id];
            if (side_[id] == 0) {
                --left;
                ++right;
            } else {
                --right;
                ++left;
            }
            const bool new_cut = left > 0 && right > 0;
            gain += db_.nets[net_id].weight *
                    (static_cast<Real>(old_cut) - static_cast<Real>(new_cut));
        }
        return gain;
    }

    std::vector<int> connectivity_order(const std::vector<int>& nodes) {
        ++stamp_;
        for (int id : nodes) membership_[id] = stamp_;
        std::vector<unsigned char> visited(db_.nodes.size(), 0);
        std::vector<int> order;
        order.reserve(nodes.size());
        std::vector<int> stack;
        for (int seed : nodes) {
            if (visited[seed]) continue;
            visited[seed] = 1;
            stack.push_back(seed);
            while (!stack.empty()) {
                const int id = stack.back();
                stack.pop_back();
                order.push_back(id);
                for (auto net_it = node_nets_[id].rbegin();
                     net_it != node_nets_[id].rend(); ++net_it) {
                    const Net& net = db_.nets[*net_it];
                    if (net.pin_count > static_cast<std::size_t>(config_.degree_limit)) {
                        continue;
                    }
                    for (std::size_t offset = net.pin_count; offset > 0; --offset) {
                        const int neighbor =
                            db_.pins[net.pin_begin + offset - 1].node;
                        if (membership_[neighbor] == stamp_ && !visited[neighbor]) {
                            visited[neighbor] = 1;
                            stack.push_back(neighbor);
                        }
                    }
                }
            }
        }
        return order;
    }

    Real node_hpwl_delta(int node_id, Real new_x, Real new_y) const {
        Real delta = 0.0;
        for (int net_id : node_nets_[node_id]) {
            const Net& net = db_.nets[net_id];
            if (net.pin_count < 2) continue;
            Real old_min_x = std::numeric_limits<Real>::infinity();
            Real old_max_x = -std::numeric_limits<Real>::infinity();
            Real old_min_y = std::numeric_limits<Real>::infinity();
            Real old_max_y = -std::numeric_limits<Real>::infinity();
            Real new_min_x = old_min_x;
            Real new_max_x = old_max_x;
            Real new_min_y = old_min_y;
            Real new_max_y = old_max_y;
            for (std::size_t p = net.pin_begin;
                 p < net.pin_begin + net.pin_count; ++p) {
                const Pin& pin = db_.pins[p];
                const Node& node = db_.nodes[pin.node];
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
            delta += net.weight *
                ((new_max_x - new_min_x) + (new_max_y - new_min_y) -
                 (old_max_x - old_min_x) - (old_max_y - old_min_y));
        }
        return delta;
    }

    std::vector<std::pair<int, Real>> rectangle_areas(
        const Node& node, Real x, Real y) const {
        std::vector<std::pair<int, Real>> areas;
        const Real left = x - 0.5 * node.width;
        const Real right = x + 0.5 * node.width;
        const Real bottom = y - 0.5 * node.height;
        const Real top = y + 0.5 * node.height;
        const int x0 = std::clamp(static_cast<int>(std::floor(
            (left - db_.xl) / density_.bin_width())),
            0, density_.bins_x() - 1);
        const int x1 = std::clamp(static_cast<int>(std::floor(
            (std::nextafter(right, left) - db_.xl) /
                density_.bin_width())),
            0, density_.bins_x() - 1);
        const int y0 = std::clamp(static_cast<int>(std::floor(
            (bottom - db_.yl) / density_.bin_height())),
            0, density_.bins_y() - 1);
        const int y1 = std::clamp(static_cast<int>(std::floor(
            (std::nextafter(top, bottom) - db_.yl) /
                density_.bin_height())),
            0, density_.bins_y() - 1);
        for (int by = y0; by <= y1; ++by) {
            const Real bin_bottom = db_.yl + by * density_.bin_height();
            const Real overlap_y = std::max<Real>(
                0.0, std::min(top, bin_bottom + density_.bin_height()) -
                         std::max(bottom, bin_bottom));
            for (int bx = x0; bx <= x1; ++bx) {
                const Real bin_left = db_.xl + bx * density_.bin_width();
                const Real overlap_x = std::max<Real>(
                    0.0, std::min(right, bin_left + density_.bin_width()) -
                             std::max(left, bin_left));
                const Real area = overlap_x * overlap_y;
                if (area > 0.0) {
                    areas.emplace_back(by * density_.bins_x() + bx, area);
                }
            }
        }
        return areas;
    }

    void place_leaf(const Region& leaf) {
        if (leaf.nodes.empty()) return;
        std::vector<int> bins;
        for (int by = leaf.y0; by < leaf.y1; ++by) {
            for (int bx = leaf.x0; bx < leaf.x1; ++bx) {
                const int bin = by * density_.bins_x() + bx;
                if (bin_capacity(density_, bin) > 0.0) bins.push_back(bin);
            }
        }
        if (bins.empty()) {
            for (int by = leaf.y0; by < leaf.y1; ++by) {
                for (int bx = leaf.x0; bx < leaf.x1; ++bx) {
                    bins.push_back(by * density_.bins_x() + bx);
                }
            }
        }
        const bool hilbert_grid = density_.bins_x() == density_.bins_y() &&
            (density_.bins_x() & (density_.bins_x() - 1)) == 0;
        std::sort(bins.begin(), bins.end(), [&](int a, int b) {
            if (!hilbert_grid) return a < b;
            return hilbert_index(density_.bins_x(), a % density_.bins_x(),
                                 a / density_.bins_x()) <
                   hilbert_index(density_.bins_x(), b % density_.bins_x(),
                                 b / density_.bins_x());
        });
        std::vector<Real> cumulative(bins.size());
        Real total_capacity = 0.0;
        for (std::size_t i = 0; i < bins.size(); ++i) {
            total_capacity += std::max<Real>(bin_capacity(density_, bins[i]), 1.0e-12);
            cumulative[i] = total_capacity;
        }
        if (config_.leaf_hpwl_guided) {
            std::vector<int> nodes = leaf.nodes;
            std::sort(nodes.begin(), nodes.end(), [&](int lhs, int rhs) {
                if (node_nets_[lhs].size() != node_nets_[rhs].size()) {
                    return node_nets_[lhs].size() > node_nets_[rhs].size();
                }
                if (db_.nodes[lhs].area() != db_.nodes[rhs].area()) {
                    return db_.nodes[lhs].area() > db_.nodes[rhs].area();
                }
                return lhs < rhs;
            });
            for (int id : nodes) {
                Node& node = db_.nodes[id];
                const int current_bx = std::clamp(
                    static_cast<int>((node.x - db_.xl) / density_.bin_width()),
                    leaf.x0, leaf.x1 - 1);
                const int current_by = std::clamp(
                    static_cast<int>((node.y - db_.yl) / density_.bin_height()),
                    leaf.y0, leaf.y1 - 1);
                int best_bin = -1;
                Real best_x = node.x;
                Real best_y = node.y;
                Real best_violation = std::numeric_limits<Real>::infinity();
                Real best_hpwl = std::numeric_limits<Real>::infinity();
                Real best_distance = std::numeric_limits<Real>::infinity();
                std::vector<std::pair<int, Real>> best_areas;
                for (int bin : bins) {
                    const int bx = bin % density_.bins_x();
                    const int by = bin / density_.bins_x();
                    if (config_.leaf_local_radius >= 0 &&
                        (std::abs(bx - current_bx) > config_.leaf_local_radius ||
                         std::abs(by - current_by) > config_.leaf_local_radius)) {
                        continue;
                    }
                    const Real x = std::clamp(
                        db_.xl + (bx + 0.5) * density_.bin_width(),
                        db_.xl + 0.5 * node.width,
                        db_.xh - 0.5 * node.width);
                    const Real y = std::clamp(
                        db_.yl + (by + 0.5) * density_.bin_height(),
                        db_.yl + 0.5 * node.height,
                        db_.yh - 0.5 * node.height);
                    const auto areas = rectangle_areas(node, x, y);
                    Real violation = 0.0;
                    for (const auto& [occupied_bin, area] : areas) {
                        violation += std::max<Real>(
                            area - remaining_capacity_[occupied_bin], 0.0);
                    }
                    const Real hpwl_delta = node_hpwl_delta(id, x, y);
                    const Real dx = (x - node.x) / density_.bin_width();
                    const Real dy = (y - node.y) / density_.bin_height();
                    const Real distance = dx * dx + dy * dy;
                    if (violation < best_violation - 1.0e-12 ||
                        (std::abs(violation - best_violation) <= 1.0e-12 &&
                         (hpwl_delta < best_hpwl - 1.0e-12 ||
                          (std::abs(hpwl_delta - best_hpwl) <= 1.0e-12 &&
                           distance < best_distance)))) {
                        best_bin = bin;
                        best_x = x;
                        best_y = y;
                        best_violation = violation;
                        best_hpwl = hpwl_delta;
                        best_distance = distance;
                        best_areas = areas;
                    }
                }
                if (best_bin < 0) continue;
                node.x = best_x;
                node.y = best_y;
                for (const auto& [bin, area] : best_areas) {
                    remaining_capacity_[bin] -= area;
                }
            }
            return;
        }
        if (config_.leaf_nearest_capacity) {
            std::vector<Real> remaining;
            remaining.reserve(bins.size());
            for (int bin : bins) remaining.push_back(bin_capacity(density_, bin));
            std::vector<int> nodes = leaf.nodes;
            std::sort(nodes.begin(), nodes.end(), [&](int lhs, int rhs) {
                const Real left = db_.nodes[lhs].area();
                const Real right = db_.nodes[rhs].area();
                if (left != right) return left > right;
                return lhs < rhs;
            });
            for (int id : nodes) {
                Node& node = db_.nodes[id];
                int best = -1;
                Real best_distance = std::numeric_limits<Real>::infinity();
                for (int pass = 0; pass < 2 && best < 0; ++pass) {
                    Real best_remaining = -std::numeric_limits<Real>::infinity();
                    for (int i = 0; i < static_cast<int>(bins.size()); ++i) {
                        const bool fits = remaining[i] + 1.0e-12 >= node.area();
                        if (pass == 0 && !fits) continue;
                        const int bx = bins[i] % density_.bins_x();
                        const int by = bins[i] / density_.bins_x();
                        const Real x = db_.xl + (bx + 0.5) * density_.bin_width();
                        const Real y = db_.yl + (by + 0.5) * density_.bin_height();
                        const Real dx = (x - node.x) / density_.bin_width();
                        const Real dy = (y - node.y) / density_.bin_height();
                        const Real distance = dx * dx + dy * dy;
                        if ((pass == 0 &&
                             (distance < best_distance ||
                              (distance == best_distance &&
                               (best < 0 || bins[i] < bins[best])))) ||
                            (pass == 1 &&
                             (remaining[i] > best_remaining ||
                              (remaining[i] == best_remaining &&
                               (distance < best_distance ||
                                (distance == best_distance &&
                                 (best < 0 || bins[i] < bins[best]))))))) {
                            best = i;
                            best_distance = distance;
                            best_remaining = remaining[i];
                        }
                    }
                }
                const int bin = bins[best];
                const int bx = bin % density_.bins_x();
                const int by = bin / density_.bins_x();
                node.x = std::clamp(
                    db_.xl + (bx + 0.5) * density_.bin_width(),
                    db_.xl + 0.5 * node.width, db_.xh - 0.5 * node.width);
                node.y = std::clamp(
                    db_.yl + (by + 0.5) * density_.bin_height(),
                    db_.yl + 0.5 * node.height, db_.yh - 0.5 * node.height);
                remaining[best] -= node.area();
            }
            return;
        }
        std::vector<int> order;
        if (config_.position_seeded) {
            order = leaf.nodes;
            std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
                const int lhs_x = std::clamp(
                    static_cast<int>((db_.nodes[lhs].x - db_.xl) /
                                     density_.bin_width()),
                    leaf.x0, leaf.x1 - 1);
                const int lhs_y = std::clamp(
                    static_cast<int>((db_.nodes[lhs].y - db_.yl) /
                                     density_.bin_height()),
                    leaf.y0, leaf.y1 - 1);
                const int rhs_x = std::clamp(
                    static_cast<int>((db_.nodes[rhs].x - db_.xl) /
                                     density_.bin_width()),
                    leaf.x0, leaf.x1 - 1);
                const int rhs_y = std::clamp(
                    static_cast<int>((db_.nodes[rhs].y - db_.yl) /
                                     density_.bin_height()),
                    leaf.y0, leaf.y1 - 1);
                if (config_.leaf_curve_order && hilbert_grid) {
                    const auto lhs_key = hilbert_index(
                        density_.bins_x(), lhs_x, lhs_y);
                    const auto rhs_key = hilbert_index(
                        density_.bins_x(), rhs_x, rhs_y);
                    if (lhs_key != rhs_key) return lhs_key < rhs_key;
                    return lhs < rhs;
                }
                const auto lhs_key = std::pair<int, int>{lhs_y, lhs_x};
                const auto rhs_key = std::pair<int, int>{rhs_y, rhs_x};
                if (lhs_key != rhs_key) return lhs_key < rhs_key;
                return lhs < rhs;
            });
        } else {
            order = connectivity_order(leaf.nodes);
        }
        Real total_area = 0.0;
        for (int id : order) total_area += db_.nodes[id].area();
        Real prefix = 0.0;
        for (int id : order) {
            Node& node = db_.nodes[id];
            const Real quantile = total_area > 0.0
                ? (prefix + 0.5 * node.area()) / total_area * total_capacity
                : 0.0;
            const auto found = std::lower_bound(cumulative.begin(), cumulative.end(), quantile);
            const std::size_t position = std::min<std::size_t>(
                std::distance(cumulative.begin(), found), bins.size() - 1);
            const int bin = bins[position];
            const int bx = bin % density_.bins_x();
            const int by = bin / density_.bins_x();
            node.x = std::clamp(
                db_.xl + (bx + 0.5) * density_.bin_width(),
                db_.xl + 0.5 * node.width, db_.xh - 0.5 * node.width);
            node.y = std::clamp(
                db_.yl + (by + 0.5) * density_.bin_height(),
                db_.yl + 0.5 * node.height, db_.yh - 0.5 * node.height);
            prefix += node.area();
        }
    }

    void measure_cut_weight() {
        for (const Net& net : db_.nets) {
            int first_leaf = -1;
            bool cut = false;
            for (std::size_t p = net.pin_begin;
                 p < net.pin_begin + net.pin_count; ++p) {
                const int node = db_.pins[p].node;
                if (db_.nodes[node].fixed) continue;
                const int leaf = leaf_of_[node];
                if (first_leaf < 0) first_leaf = leaf;
                else if (leaf != first_leaf) {
                    cut = true;
                    break;
                }
            }
            if (cut) stats_.cut_net_weight += net.weight;
        }
    }
};

}  // namespace

BisectionStats recursive_hypergraph_bisection(
    Database& db, ExactOverlapDensity& density, const BisectionConfig& config) {
    if (config.leaf_bins <= 0 || config.degree_limit < 2 ||
        config.fm_passes < 0 || config.balance_tolerance < 0.0 ||
        config.balance_tolerance >= 0.5 || config.atomic_max_nodes < 2 ||
        config.atomic_degree_limit < 2 || config.atomic_rounds < 0 ||
        config.atomic_release_depth < 0 || config.leaf_local_radius < -1 ||
        config.max_depth < -1) {
        throw std::invalid_argument("invalid recursive bisection configuration");
    }
    BisectionStats stats;
    if (!config.enabled) return stats;
    const auto started = std::chrono::steady_clock::now();
    ExactHpwl hpwl(db);
    stats.initial_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    stats.initial_overflow = density.evaluate(0.0, 1.0, nullptr, nullptr).overflow;
    Bisector bisector(db, density, config, stats);
    bisector.run();
    stats.final_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    stats.final_overflow = density.evaluate(0.0, 1.0, nullptr, nullptr).overflow;
    stats.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return stats;
}

}  // namespace ea
