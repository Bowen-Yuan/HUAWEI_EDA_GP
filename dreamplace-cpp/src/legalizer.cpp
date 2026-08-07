#include "legalizer.h"

#include "bookshelf.h"
#include "wirelength.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <filesystem>

namespace dpcpp {
namespace {

struct Interval {
    Real lo = 0.0;
    Real hi = 0.0;
};

struct Segment {
    int row = -1;
    Real lo = 0.0;
    Real hi = 0.0;
    std::vector<Interval> free;
    std::vector<int> cells;
};

struct Desired {
    Real x = 0.0;
    Real y = 0.0;
};

Real snap_up(Real x, const Row& row) {
    return row.origin + std::ceil((x - row.origin) / row.site_spacing - 1.0e-9) *
                            row.site_spacing;
}

Real snap_nearest(Real x, const Row& row) {
    return row.origin + std::round((x - row.origin) / row.site_spacing) *
                            row.site_spacing;
}

std::vector<Segment> build_segments(const Database& db) {
    std::vector<Segment> result;
    for (int r = 0; r < static_cast<int>(db.rows.size()); ++r) {
        const Row& row = db.rows[r];
        std::vector<Interval> blocked;
        for (int id : db.fixed_ids) {
            const Node& fixed = db.nodes[id];
            if (fixed.terminal_ni) continue;
            const Real bottom = fixed.y - 0.5 * fixed.height;
            const Real top = fixed.y + 0.5 * fixed.height;
            if (top <= row.y + 1.0e-9 || bottom >= row.y + row.height - 1.0e-9) continue;
            const Real lo = std::max(row.origin, fixed.x - 0.5 * fixed.width);
            const Real hi = std::min(row.xh(), fixed.x + 0.5 * fixed.width);
            if (hi > lo) blocked.push_back({lo, hi});
        }
        std::sort(blocked.begin(), blocked.end(), [](const Interval& a, const Interval& b) {
            return a.lo < b.lo;
        });
        std::vector<Interval> merged;
        for (const Interval& interval : blocked) {
            if (merged.empty() || interval.lo > merged.back().hi) merged.push_back(interval);
            else merged.back().hi = std::max(merged.back().hi, interval.hi);
        }
        Real cursor = row.origin;
        for (const Interval& interval : merged) {
            const Real end = std::min(row.xh(), interval.lo);
            const Real lo = snap_up(cursor, row);
            const Real hi = row.origin + std::floor((end - row.origin) / row.site_spacing + 1.0e-9) *
                                           row.site_spacing;
            if (hi > lo + 1.0e-9) result.push_back({r, lo, hi, {{lo, hi}}, {}});
            cursor = std::max(cursor, interval.hi);
        }
        const Real lo = snap_up(cursor, row);
        const Real hi = row.origin + std::floor((row.xh() - row.origin) / row.site_spacing + 1.0e-9) *
                                       row.site_spacing;
        if (hi > lo + 1.0e-9) result.push_back({r, lo, hi, {{lo, hi}}, {}});
    }
    return result;
}

void save_legal_snapshot(const Database& db, const LegalizeConfig& config,
                         const char* stage) {
    if (config.snapshot_dir.empty()) return;
    write_bookshelf_pl(db, (std::filesystem::path(config.snapshot_dir) /
                            (std::string("legal_") + stage + ".pl")).string());
}

bool fit_in_segment(const Segment& segment, const Row& row, Real target_left,
                    Real width, Real& placed_left, int& interval_index) {
    Real best_cost = std::numeric_limits<Real>::infinity();
    bool found = false;
    for (int i = 0; i < static_cast<int>(segment.free.size()); ++i) {
        const Interval& space = segment.free[i];
        if (space.hi - space.lo + 1.0e-9 < width) continue;
        Real x = std::clamp(snap_nearest(target_left, row), space.lo, space.hi - width);
        x = snap_up(x, row);
        if (x + width > space.hi + 1.0e-9) x -= row.site_spacing;
        if (x < space.lo - 1.0e-9 || x + width > space.hi + 1.0e-9) continue;
        const Real cost = std::abs(x - target_left);
        if (cost < best_cost) {
            best_cost = cost;
            placed_left = x;
            interval_index = i;
            found = true;
        }
    }
    return found;
}

void occupy(Segment& segment, int interval_index, Real lo, Real hi) {
    const Interval old = segment.free[interval_index];
    segment.free.erase(segment.free.begin() + interval_index);
    if (hi < old.hi - 1.0e-9) segment.free.insert(
        segment.free.begin() + interval_index, {hi, old.hi});
    if (old.lo < lo - 1.0e-9) segment.free.insert(
        segment.free.begin() + interval_index, {old.lo, lo});
}

int nearest_row(const Database& db, Real target_bottom) {
    auto it = std::lower_bound(db.rows.begin(), db.rows.end(), target_bottom,
        [](const Row& row, Real y) { return row.y < y; });
    int index = static_cast<int>(it - db.rows.begin());
    if (index == static_cast<int>(db.rows.size())) --index;
    if (index > 0 && std::abs(db.rows[index - 1].y - target_bottom) <
                     std::abs(db.rows[index].y - target_bottom)) --index;
    return std::max(0, index);
}

void greedy_legalize(Database& db, const std::vector<Desired>& desired,
                     std::vector<Segment>& segments, int row_search_limit) {
    std::vector<std::vector<int>> segments_by_row(db.rows.size());
    for (int s = 0; s < static_cast<int>(segments.size()); ++s) {
        segments_by_row[segments[s].row].push_back(s);
    }
    std::vector<int> order = db.movable_ids;
    // Preserve the global placement's left-to-right order.  A global
    // width-first order is robust against fragmentation, but destroys
    // locality by letting a small set of wide cells consume every nearby
    // row before the spatial sweep reaches ordinary cells.
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
        if (desired[a].x != desired[b].x) return desired[a].x < desired[b].x;
        if (db.nodes[a].width != db.nodes[b].width)
            return db.nodes[a].width > db.nodes[b].width;
        return desired[a].y < desired[b].y;
    });

    int placed = 0;
    for (int id : order) {
        Node& node = db.nodes[id];
        const int base_row = nearest_row(db, desired[id].y - 0.5 * node.height);
        Real best_cost = std::numeric_limits<Real>::infinity();
        Real best_left = 0.0;
        int best_segment = -1;
        int best_interval = -1;
        const int maximum_radius = std::max(row_search_limit,
                                             static_cast<int>(db.rows.size()));
        for (int radius = 0; radius <= maximum_radius; ++radius) {
            bool tested = false;
            for (int direction : {-1, 1}) {
                if (radius == 0 && direction == 1) continue;
                const int r = base_row + direction * radius;
                if (r < 0 || r >= static_cast<int>(db.rows.size())) continue;
                const Row& row = db.rows[r];
                if (node.height > row.height + 1.0e-6) continue;
                tested = true;
                for (int s : segments_by_row[r]) {
                    Real left = 0.0;
                    int interval = -1;
                    if (!fit_in_segment(segments[s], row,
                                        desired[id].x - 0.5 * node.width,
                                        node.width, left, interval)) continue;
                    const Real cost = std::abs(left + 0.5 * node.width - desired[id].x) +
                                      2.0 * std::abs(row.y + 0.5 * node.height - desired[id].y);
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_left = left;
                        best_segment = s;
                        best_interval = interval;
                    }
                }
            }
            if (best_segment >= 0 && (radius >= row_search_limit ||
                2.0 * radius * db.rows.front().height > best_cost)) break;
            if (!tested && radius >= static_cast<int>(db.rows.size())) break;
        }
        if (best_segment < 0) {
            throw std::runtime_error("greedy legalization cannot place node " + node.name);
        }
        Segment& segment = segments[best_segment];
        const Row& row = db.rows[segment.row];
        node.x = best_left + 0.5 * node.width;
        node.y = row.y + 0.5 * node.height;
        occupy(segment, best_interval, best_left, best_left + node.width);
        segment.cells.push_back(id);
        ++placed;
    }
    std::cout << "[Greedy] placed=" << placed << " segments=" << segments.size() << '\n';
}

void abacus_segment(Database& db, Segment& segment,
                    const std::vector<Desired>& desired) {
    if (segment.cells.empty()) return;
    const Row& row = db.rows[segment.row];
    std::sort(segment.cells.begin(), segment.cells.end(), [&](int a, int b) {
        return db.nodes[a].x < db.nodes[b].x;
    });
    const int n = static_cast<int>(segment.cells.size());
    std::vector<Real> prefix(n + 1, 0.0);
    for (int i = 0; i < n; ++i) prefix[i + 1] = prefix[i] + db.nodes[segment.cells[i]].width;
    const Real zlo = segment.lo;
    const Real zhi = segment.hi - prefix[n];
    if (zhi < zlo - 1.0e-6) throw std::runtime_error("overfull Abacus segment");

    struct Block { int begin; int end; Real sum; int count; };
    std::vector<Block> blocks;
    for (int i = 0; i < n; ++i) {
        Real target = desired[segment.cells[i]].x - 0.5 * db.nodes[segment.cells[i]].width - prefix[i];
        target = std::clamp(target, zlo, zhi);
        blocks.push_back({i, i + 1, target, 1});
        while (blocks.size() >= 2) {
            const Block& a = blocks[blocks.size() - 2];
            const Block& b = blocks.back();
            if (a.sum / a.count <= b.sum / b.count + 1.0e-12) break;
            Block merged{a.begin, b.end, a.sum + b.sum, a.count + b.count};
            blocks.pop_back();
            blocks.back() = merged;
        }
    }
    std::vector<Real> z(n);
    for (const Block& block : blocks) {
        Real value = std::clamp(block.sum / block.count, zlo, zhi);
        value = std::clamp(snap_nearest(value, row), zlo, zhi);
        for (int i = block.begin; i < block.end; ++i) z[i] = value;
    }
    for (int i = 1; i < n; ++i) z[i] = std::max(z[i], z[i - 1]);
    for (int i = n - 2; i >= 0; --i) z[i] = std::min(z[i], z[i + 1]);
    for (int i = 0; i < n; ++i) {
        Node& node = db.nodes[segment.cells[i]];
        const Real left = z[i] + prefix[i];
        node.x = left + 0.5 * node.width;
        node.y = row.y + 0.5 * node.height;
    }
}

Real affected_hpwl(const Database& db, const std::vector<int>& nets) {
    Real total = 0.0;
    for (int net_id : nets) {
        const Net& net = db.nets[net_id];
        if (net.pins.size() < 2) continue;
        Real xmin = std::numeric_limits<Real>::infinity();
        Real xmax = -xmin;
        Real ymin = xmin;
        Real ymax = -xmin;
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            xmin = std::min(xmin, x); xmax = std::max(xmax, x);
            ymin = std::min(ymin, y); ymax = std::max(ymax, y);
        }
        total += net.weight * ((xmax - xmin) + (ymax - ymin));
    }
    return total;
}

std::vector<std::vector<int>> build_incident_nets(const Database& db) {
    std::vector<std::vector<int>> incident(db.nodes.size());
    for (int e = 0; e < static_cast<int>(db.nets.size()); ++e) {
        for (const Pin& pin : db.nets[e].pins) incident[pin.node].push_back(e);
    }
    return incident;
}

std::vector<int> union_incident(const std::vector<std::vector<int>>& incident,
                                const std::vector<int>& cells) {
    std::vector<int> nets;
    for (int id : cells) nets.insert(nets.end(), incident[id].begin(), incident[id].end());
    std::sort(nets.begin(), nets.end());
    nets.erase(std::unique(nets.begin(), nets.end()), nets.end());
    return nets;
}

void detailed_reorder(Database& db, std::vector<Segment>& segments, int passes) {
    const std::vector<std::vector<int>> incident = build_incident_nets(db);
    std::size_t accepted_total = 0;
    for (int pass = 0; pass < passes; ++pass) {
        std::size_t accepted = 0;
        for (Segment& segment : segments) {
            std::sort(segment.cells.begin(), segment.cells.end(), [&](int a, int b) {
                return db.nodes[a].x < db.nodes[b].x;
            });
            for (int i = 0; i + 1 < static_cast<int>(segment.cells.size()); ++i) {
                const int a = segment.cells[i];
                const int b = segment.cells[i + 1];
                Node& na = db.nodes[a];
                Node& nb = db.nodes[b];
                const Real left = na.x - 0.5 * na.width;
                const Real gap = (nb.x - 0.5 * nb.width) - (na.x + 0.5 * na.width);
                std::vector<int> nets = incident[a];
                nets.insert(nets.end(), incident[b].begin(), incident[b].end());
                std::sort(nets.begin(), nets.end());
                nets.erase(std::unique(nets.begin(), nets.end()), nets.end());
                const Real before = affected_hpwl(db, nets);
                const Real old_ax = na.x;
                const Real old_bx = nb.x;
                nb.x = left + 0.5 * nb.width;
                na.x = left + nb.width + gap + 0.5 * na.width;
                const Real after = affected_hpwl(db, nets);
                if (after + 1.0e-9 < before) {
                    std::swap(segment.cells[i], segment.cells[i + 1]);
                    ++accepted;
                } else {
                    na.x = old_ax;
                    nb.x = old_bx;
                }
            }
        }
        accepted_total += accepted;
        std::cout << "[DP] pass=" << pass << " accepted_swaps=" << accepted << '\n';
        if (accepted == 0) break;
    }
    std::cout << "[DP] accepted_total=" << accepted_total << '\n';
}

void assign_window(Database& db, const std::vector<int>& order, Real left,
                   const std::vector<Real>& gaps) {
    Real cursor = left;
    for (int i = 0; i < static_cast<int>(order.size()); ++i) {
        Node& node = db.nodes[order[i]];
        node.x = cursor + 0.5 * node.width;
        cursor += node.width;
        if (i < static_cast<int>(gaps.size())) cursor += gaps[i];
    }
}

std::size_t k_reorder(Database& db, std::vector<Segment>& segments,
                      int window_size, int passes) {
    const std::vector<std::vector<int>> incident = build_incident_nets(db);
    std::size_t accepted_total = 0;
    for (int pass = 0; pass < passes; ++pass) {
        std::size_t accepted = 0;
        for (Segment& segment : segments) {
            std::sort(segment.cells.begin(), segment.cells.end(), [&](int a, int b) {
                return db.nodes[a].x < db.nodes[b].x;
            });
            const int count = static_cast<int>(segment.cells.size());
            for (int begin = 0; begin + 1 < count;
                 begin += std::max(1, window_size - 1)) {
                const int size = std::min(window_size, count - begin);
                if (size < 2) continue;
                std::vector<int> original(segment.cells.begin() + begin,
                                          segment.cells.begin() + begin + size);
                const Real left = db.nodes[original.front()].x -
                                  0.5 * db.nodes[original.front()].width;
                std::vector<Real> gaps(size - 1, 0.0);
                for (int i = 0; i + 1 < size; ++i) {
                    gaps[i] = (db.nodes[original[i + 1]].x -
                               0.5 * db.nodes[original[i + 1]].width) -
                              (db.nodes[original[i]].x +
                               0.5 * db.nodes[original[i]].width);
                }
                const std::vector<int> nets = union_incident(incident, original);
                const Real before = affected_hpwl(db, nets);
                Real best = before;
                std::vector<int> best_order = original;
                std::vector<int> candidate = original;
                std::sort(candidate.begin(), candidate.end());
                do {
                    if (candidate == original) continue;
                    assign_window(db, candidate, left, gaps);
                    const Real value = affected_hpwl(db, nets);
                    if (value + 1.0e-9 < best) {
                        best = value;
                        best_order = candidate;
                    }
                } while (std::next_permutation(candidate.begin(), candidate.end()));
                assign_window(db, best_order, left, gaps);
                if (best + 1.0e-9 < before) {
                    std::copy(best_order.begin(), best_order.end(),
                              segment.cells.begin() + begin);
                    ++accepted;
                } else {
                    assign_window(db, original, left, gaps);
                }
            }
        }
        accepted_total += accepted;
        std::cout << "[K-Reorder] pass=" << pass << " accepted=" << accepted << '\n';
        if (accepted == 0) break;
    }
    return accepted_total;
}

std::pair<Real, Real> net_target(const Database& db, int cell,
                                 const std::vector<int>& nets) {
    Real tx = 0.0, ty = 0.0, total_weight = 0.0;
    for (int net_id : nets) {
        const Net& net = db.nets[net_id];
        Real xmin = std::numeric_limits<Real>::infinity();
        Real xmax = -xmin, ymin = xmin, ymax = -xmin;
        bool found = false;
        for (const Pin& pin : net.pins) {
            if (pin.node == cell) continue;
            const Node& node = db.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            xmin = std::min(xmin, x); xmax = std::max(xmax, x);
            ymin = std::min(ymin, y); ymax = std::max(ymax, y);
            found = true;
        }
        if (found) {
            tx += net.weight * 0.5 * (xmin + xmax);
            ty += net.weight * 0.5 * (ymin + ymax);
            total_weight += net.weight;
        }
    }
    if (total_weight == 0.0) return {db.nodes[cell].x, db.nodes[cell].y};
    return {tx / total_weight, ty / total_weight};
}

using SizeKey = std::pair<long long, long long>;

SizeKey size_key(const Node& node) {
    constexpr Real scale = 1.0e6;
    return {std::llround(node.width * scale), std::llround(node.height * scale)};
}

std::map<SizeKey, std::vector<int>> equal_size_groups(const Database& db) {
    std::map<SizeKey, std::vector<int>> groups;
    for (int id : db.movable_ids) groups[size_key(db.nodes[id])].push_back(id);
    for (auto& entry : groups) {
        std::sort(entry.second.begin(), entry.second.end(), [&](int a, int b) {
            return db.nodes[a].x < db.nodes[b].x;
        });
    }
    return groups;
}

std::size_t global_swap(Database& db, int passes) {
    const std::vector<std::vector<int>> incident = build_incident_nets(db);
    std::size_t accepted_total = 0;
    for (int pass = 0; pass < passes; ++pass) {
        auto groups = equal_size_groups(db);
        std::unordered_set<int> used;
        std::size_t accepted = 0;
        for (auto& entry : groups) {
            std::vector<int>& group = entry.second;
            if (group.size() < 2) continue;
            for (int slot_index = 0; slot_index < static_cast<int>(group.size());
                 ++slot_index) {
                const int id = group[slot_index];
                if (used.count(id)) continue;
                const auto target = net_target(db, id, incident[id]);
                auto position = std::lower_bound(group.begin(), group.end(), target.first,
                    [&](int candidate, Real x) { return db.nodes[candidate].x < x; });
                const int center = static_cast<int>(position - group.begin());
                int best_id = -1;
                int best_index = -1;
                Real best_distance = std::numeric_limits<Real>::infinity();
                for (int j = std::max(0, center - 8);
                     j < std::min(static_cast<int>(group.size()), center + 9); ++j) {
                    const int candidate = group[j];
                    if (candidate == id || used.count(candidate)) continue;
                    const Real distance = std::abs(db.nodes[candidate].x - target.first) +
                                          std::abs(db.nodes[candidate].y - target.second);
                    if (distance < best_distance) {
                        best_distance = distance;
                        best_id = candidate;
                        best_index = j;
                    }
                }
                if (best_id < 0) continue;
                const std::vector<int> cells{id, best_id};
                const std::vector<int> nets = union_incident(incident, cells);
                const Real before = affected_hpwl(db, nets);
                std::swap(db.nodes[id].x, db.nodes[best_id].x);
                std::swap(db.nodes[id].y, db.nodes[best_id].y);
                const Real after = affected_hpwl(db, nets);
                if (after + 1.0e-9 < before) {
                    std::swap(group[slot_index], group[best_index]);
                    used.insert(id); used.insert(best_id); ++accepted;
                } else {
                    std::swap(db.nodes[id].x, db.nodes[best_id].x);
                    std::swap(db.nodes[id].y, db.nodes[best_id].y);
                }
            }
        }
        accepted_total += accepted;
        std::cout << "[GlobalSwap] pass=" << pass << " accepted=" << accepted << '\n';
        if (accepted == 0) break;
    }
    return accepted_total;
}

bool disjoint_nets(const std::vector<int>& nets,
                   const std::unordered_set<int>& occupied) {
    for (int net : nets) if (occupied.count(net)) return false;
    return true;
}

std::size_t independent_set_matching(Database& db, int maximum_size) {
    const std::vector<std::vector<int>> incident = build_incident_nets(db);
    auto groups = equal_size_groups(db);
    std::size_t accepted = 0;
    for (auto& entry : groups) {
        const std::vector<int>& group = entry.second;
        if (group.size() < 2) continue;
        std::vector<char> consumed(group.size(), false);
        for (int seed_index = 0; seed_index < static_cast<int>(group.size()); ++seed_index) {
            if (consumed[seed_index]) continue;
            std::vector<int> cells{group[seed_index]};
            std::vector<int> indices{seed_index};
            std::unordered_set<int> net_set(incident[group[seed_index]].begin(),
                                            incident[group[seed_index]].end());
            const auto target = net_target(db, group[seed_index], incident[group[seed_index]]);
            auto middle = std::lower_bound(group.begin(), group.end(), target.first,
                [&](int candidate, Real x) { return db.nodes[candidate].x < x; });
            const int center = static_cast<int>(middle - group.begin());
            for (int radius = 0; radius < static_cast<int>(group.size()) &&
                                 static_cast<int>(cells.size()) < maximum_size; ++radius) {
                for (int sign : {-1, 1}) {
                    if (radius == 0 && sign == 1) continue;
                    const int j = center + sign * radius;
                    if (j < 0 || j >= static_cast<int>(group.size()) || consumed[j] ||
                        j == seed_index ||
                        !disjoint_nets(incident[group[j]], net_set)) continue;
                    cells.push_back(group[j]); indices.push_back(j);
                    net_set.insert(incident[group[j]].begin(), incident[group[j]].end());
                    if (static_cast<int>(cells.size()) == maximum_size) break;
                }
            }
            if (cells.size() < 2) { consumed[seed_index] = true; continue; }
            for (int index : indices) consumed[index] = true;
            std::vector<std::pair<Real, Real>> slots;
            for (int id : cells) slots.push_back({db.nodes[id].x, db.nodes[id].y});
            std::vector<int> permutation(cells.size());
            std::iota(permutation.begin(), permutation.end(), 0);
            Real best_cost = std::numeric_limits<Real>::infinity();
            std::vector<int> best_permutation = permutation;
            do {
                Real cost = 0.0;
                for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
                    Node& node = db.nodes[cells[i]];
                    const Real old_x = node.x, old_y = node.y;
                    node.x = slots[permutation[i]].first;
                    node.y = slots[permutation[i]].second;
                    cost += affected_hpwl(db, incident[cells[i]]);
                    node.x = old_x; node.y = old_y;
                }
                if (cost < best_cost) { best_cost = cost; best_permutation = permutation; }
            } while (std::next_permutation(permutation.begin(), permutation.end()));
            const std::vector<int> nets = union_incident(incident, cells);
            const Real before = affected_hpwl(db, nets);
            for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
                db.nodes[cells[i]].x = slots[best_permutation[i]].first;
                db.nodes[cells[i]].y = slots[best_permutation[i]].second;
            }
            const Real after = affected_hpwl(db, nets);
            if (after + 1.0e-9 < before) {
                ++accepted;
            } else {
                for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
                    db.nodes[cells[i]].x = slots[i].first;
                    db.nodes[cells[i]].y = slots[i].second;
                }
            }
        }
    }
    std::cout << "[IndependentSet] accepted=" << accepted << '\n';
    return accepted;
}

bool overlaps(Real alo, Real ahi, Real blo, Real bhi, Real tolerance) {
    return std::min(ahi, bhi) - std::max(alo, blo) > tolerance;
}

}  // namespace

LegalityResult check_legality(const Database& db, Real tolerance) {
    LegalityResult result;
    std::vector<std::vector<int>> by_row(db.rows.size());
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real left = node.x - 0.5 * node.width;
        const Real right = node.x + 0.5 * node.width;
        const Real bottom = node.y - 0.5 * node.height;
        int row_index = nearest_row(db, bottom);
        const Row& row = db.rows[row_index];
        bool boundary_ok = left >= row.origin - tolerance && right <= row.xh() + tolerance &&
                           bottom >= row.y - tolerance &&
                           bottom + node.height <= row.y + row.height + tolerance;
        if (!boundary_ok) {
            ++result.boundary_errors;
            if (result.first_error.empty()) result.first_error = "boundary: " + node.name;
        }
        const Real site = (left - row.origin) / row.site_spacing;
        if (std::abs(site - std::round(site)) > tolerance) {
            ++result.alignment_errors;
            if (result.first_error.empty()) result.first_error = "site alignment: " + node.name;
        }
        by_row[row_index].push_back(id);
        for (int fixed_id : db.fixed_ids) {
            const Node& fixed = db.nodes[fixed_id];
            if (fixed.terminal_ni) continue;
            if (overlaps(left, right, fixed.x - 0.5 * fixed.width,
                         fixed.x + 0.5 * fixed.width, tolerance) &&
                overlaps(bottom, bottom + node.height, fixed.y - 0.5 * fixed.height,
                         fixed.y + 0.5 * fixed.height, tolerance)) {
                ++result.overlap_errors;
                if (result.first_error.empty()) result.first_error =
                    "fixed overlap: " + node.name + " / " + fixed.name;
                break;
            }
        }
    }
    for (std::vector<int>& ids : by_row) {
        std::sort(ids.begin(), ids.end(), [&](int a, int b) {
            return db.nodes[a].x - 0.5 * db.nodes[a].width <
                   db.nodes[b].x - 0.5 * db.nodes[b].width;
        });
        for (int i = 0; i + 1 < static_cast<int>(ids.size()); ++i) {
            const Node& a = db.nodes[ids[i]];
            const Node& b = db.nodes[ids[i + 1]];
            if (a.x + 0.5 * a.width > b.x - 0.5 * b.width + tolerance) {
                ++result.overlap_errors;
                if (result.first_error.empty()) result.first_error =
                    "movable overlap: " + a.name + " / " + b.name;
            }
        }
    }
    result.legal = result.boundary_errors == 0 && result.alignment_errors == 0 &&
                   result.overlap_errors == 0;
    return result;
}

LegalizeResult legalize_and_refine(Database& db, const LegalizeConfig& config) {
    LegalizeResult result;
    result.hpwl_before = exact_hpwl(db);
    std::vector<Desired> desired(db.nodes.size());
    for (int id : db.movable_ids) desired[id] = {db.nodes[id].x, db.nodes[id].y};
    std::vector<Segment> segments = build_segments(db);
    greedy_legalize(db, desired, segments, config.row_search_limit);
    result.hpwl_after_greedy = exact_hpwl(db);
    save_legal_snapshot(db, config, "greedy");
    LegalityResult legal = check_legality(db);
    if (!legal.legal) throw std::runtime_error("Greedy produced illegal placement: " + legal.first_error);

    if (config.run_abacus) {
        for (Segment& segment : segments) abacus_segment(db, segment, desired);
    }
    result.hpwl_after_abacus = exact_hpwl(db);
    save_legal_snapshot(db, config, "abacus");
    legal = check_legality(db);
    if (!legal.legal) throw std::runtime_error("Abacus produced illegal placement: " + legal.first_error);

    if (config.run_detailed) detailed_reorder(db, segments, config.detailed_passes);
    save_legal_snapshot(db, config, "adjacent");
    if (config.run_dreamplace_detailed) {
        k_reorder(db, segments, config.k_reorder_size, config.detailed_passes);
    }
    result.hpwl_after_k_reorder = exact_hpwl(db);
    save_legal_snapshot(db, config, "k_reorder");
    legal = check_legality(db);
    if (!legal.legal) throw std::runtime_error(
        "K-Reorder produced illegal placement: " + legal.first_error);

    if (config.run_dreamplace_detailed) {
        global_swap(db, config.global_swap_passes);
    }
    result.hpwl_after_global_swap = exact_hpwl(db);
    save_legal_snapshot(db, config, "global_swap");
    legal = check_legality(db);
    if (!legal.legal) throw std::runtime_error(
        "global swap produced illegal placement: " + legal.first_error);

    if (config.run_dreamplace_detailed) {
        independent_set_matching(db, config.independent_set_size);
    }
    result.hpwl_after_independent_set = exact_hpwl(db);
    save_legal_snapshot(db, config, "independent_set");
    result.hpwl_after_detailed = result.hpwl_after_independent_set;
    result.legality = check_legality(db);
    if (!result.legality.legal) {
        throw std::runtime_error("detailed placement is illegal: " + result.legality.first_error);
    }
    std::cout << "[Legalize] HPWL before=" << result.hpwl_before
              << " greedy=" << result.hpwl_after_greedy
              << " abacus=" << result.hpwl_after_abacus
              << " k_reorder=" << result.hpwl_after_k_reorder
              << " global_swap=" << result.hpwl_after_global_swap
              << " independent_set=" << result.hpwl_after_independent_set
              << " detailed=" << result.hpwl_after_detailed << '\n';
    return result;
}

}  // namespace dpcpp
