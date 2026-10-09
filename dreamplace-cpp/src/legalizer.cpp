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

std::vector<std::vector<int>> index_fixed_obstacles_by_row(const Database& db) {
    std::vector<std::vector<int>> result(db.rows.size());
    if (db.rows.empty()) return result;
    Real maximum_row_height = 0.0;
    for (const Row& row : db.rows)
        maximum_row_height = std::max(maximum_row_height, row.height);
    for (int id : db.fixed_ids) {
        const Node& fixed = db.nodes[id];
        if (fixed.terminal_ni) continue;
        const Real bottom = fixed.y - 0.5 * fixed.height;
        const Real top = fixed.y + 0.5 * fixed.height;
        auto first = std::lower_bound(
            db.rows.begin(), db.rows.end(), bottom - maximum_row_height,
            [](const Row& row, Real y) { return row.y < y; });
        for (auto it = first; it != db.rows.end() && it->y < top - 1.0e-9; ++it) {
            if (it->y + it->height <= bottom + 1.0e-9) continue;
            result[static_cast<std::size_t>(it - db.rows.begin())].push_back(id);
        }
    }
    return result;
}

std::vector<Segment> build_segments(const Database& db) {
    std::vector<Segment> result;
    const std::vector<std::vector<int>> obstacles_by_row =
        index_fixed_obstacles_by_row(db);
    for (int r = 0; r < static_cast<int>(db.rows.size()); ++r) {
        const Row& row = db.rows[r];
        std::vector<Interval> blocked;
        for (int id : obstacles_by_row[r]) {
            const Node& fixed = db.nodes[id];
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

int covered_row_count(const Database& db, int base_row, const Node& node) {
    if (base_row < 0 || base_row >= static_cast<int>(db.rows.size())) return 0;
    const Real bottom = db.rows[base_row].y;
    const Real top = bottom + node.height;
    int count = 0;
    Real cursor = bottom;
    for (int r = base_row; r < static_cast<int>(db.rows.size()) &&
         cursor < top - 1.0e-6; ++r) {
        const Row& row = db.rows[r];
        if (std::abs(row.y - cursor) > 1.0e-6) return 0;
        cursor = row.y + row.height;
        ++count;
    }
    return std::abs(cursor - top) <= 1.0e-6 ? count : 0;
}

bool fit_multiline_cell(const Database& db,
                        const std::vector<std::vector<int>>& segments_by_row,
                        const std::vector<Segment>& segments,
                        int base_row, const Node& node, Real target_left,
                        Real& placed_left,
                        std::vector<std::pair<int, int>>& reservations) {
    const int row_count = covered_row_count(db, base_row, node);
    if (row_count <= 1) return false;
    std::vector<Interval> common;
    for (int s : segments_by_row[base_row])
        common.insert(common.end(), segments[s].free.begin(), segments[s].free.end());
    for (int offset = 1; offset < row_count && !common.empty(); ++offset) {
        std::vector<Interval> next;
        for (const Interval& a : common) {
            for (int s : segments_by_row[base_row + offset]) {
                for (const Interval& b : segments[s].free) {
                    const Real lo = std::max(a.lo, b.lo);
                    const Real hi = std::min(a.hi, b.hi);
                    if (hi - lo + 1.0e-9 >= node.width) next.push_back({lo, hi});
                }
            }
        }
        common.swap(next);
    }
    const Row& base = db.rows[base_row];
    Real best_cost = std::numeric_limits<Real>::infinity();
    bool found = false;
    for (const Interval& space : common) {
        Real left = std::clamp(snap_nearest(target_left, base),
                               space.lo, space.hi - node.width);
        left = snap_up(left, base);
        if (left + node.width > space.hi + 1.0e-9) left -= base.site_spacing;
        if (left < space.lo - 1.0e-9 ||
            left + node.width > space.hi + 1.0e-9) continue;
        bool aligned = true;
        for (int offset = 1; offset < row_count; ++offset) {
            const Row& row = db.rows[base_row + offset];
            const Real site = (left - row.origin) / row.site_spacing;
            if (std::abs(site - std::round(site)) > 1.0e-6) {
                aligned = false;
                break;
            }
        }
        if (!aligned) continue;
        const Real cost = std::abs(left - target_left);
        if (cost < best_cost) {
            best_cost = cost;
            placed_left = left;
            found = true;
        }
    }
    if (!found) return false;

    reservations.clear();
    for (int offset = 0; offset < row_count; ++offset) {
        bool reserved = false;
        for (int s : segments_by_row[base_row + offset]) {
            for (int interval = 0;
                 interval < static_cast<int>(segments[s].free.size()); ++interval) {
                const Interval& space = segments[s].free[interval];
                if (placed_left >= space.lo - 1.0e-9 &&
                    placed_left + node.width <= space.hi + 1.0e-9) {
                    reservations.emplace_back(s, interval);
                    reserved = true;
                    break;
                }
            }
            if (reserved) break;
        }
        if (!reserved) return false;
    }
    return true;
}

std::vector<Segment> rebuild_segment_cells(const Database& db) {
    std::vector<Segment> segments = build_segments(db);
    std::vector<std::vector<int>> by_row(db.rows.size());
    for (int s = 0; s < static_cast<int>(segments.size()); ++s)
        by_row[segments[s].row].push_back(s);
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const int row_index = nearest_row(db, node.y - 0.5 * node.height);
        const Real left = node.x - 0.5 * node.width;
        const Real right = node.x + 0.5 * node.width;
        int owner = -1;
        for (int s : by_row[row_index]) {
            if (left >= segments[s].lo - 1.0e-6 &&
                right <= segments[s].hi + 1.0e-6) {
                owner = s;
                break;
            }
        }
        if (owner < 0)
            throw std::runtime_error("legal cell has no obstacle-free segment: " + node.name);
        segments[owner].cells.push_back(id);
    }
    return segments;
}

using MovableCoordinates = std::vector<std::pair<Real, Real>>;

MovableCoordinates save_movable_coordinates(const Database& db) {
    MovableCoordinates coordinates;
    coordinates.reserve(db.movable_ids.size());
    for (int id : db.movable_ids)
        coordinates.emplace_back(db.nodes[id].x, db.nodes[id].y);
    return coordinates;
}

void restore_movable_coordinates(Database& db,
                                 const MovableCoordinates& coordinates) {
    for (int i = 0; i < static_cast<int>(db.movable_ids.size()); ++i) {
        Node& node = db.nodes[db.movable_ids[i]];
        node.x = coordinates[i].first;
        node.y = coordinates[i].second;
    }
}

void greedy_legalize(Database& db, const std::vector<Desired>& desired,
                     std::vector<Segment>& segments, int row_search_limit) {
    const MovableCoordinates original = save_movable_coordinates(db);

    auto try_order = [&](std::vector<int> order) -> bool {
        // Every retry starts from the same continuous GP point and a fresh
        // obstacle-aware segment map.  This avoids carrying fragmentation
        // created by a failed greedy ordering into the fallback ordering.
        restore_movable_coordinates(db, original);
        segments = build_segments(db);
        std::vector<std::vector<int>> segments_by_row(db.rows.size());
        for (int s = 0; s < static_cast<int>(segments.size()); ++s)
            segments_by_row[segments[s].row].push_back(s);

        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            const bool tall_a = db.nodes[a].height > db.rows.front().height + 1.0e-6;
            const bool tall_b = db.nodes[b].height > db.rows.front().height + 1.0e-6;
            if (tall_a != tall_b) return tall_a;
            if (tall_a && db.nodes[a].width != db.nodes[b].width)
                return db.nodes[a].width > db.nodes[b].width;
            return false;
        });

        int placed = 0;
        for (int id : order) {
        Node& node = db.nodes[id];
        const int base_row = nearest_row(db, desired[id].y - 0.5 * node.height);
        Real best_cost = std::numeric_limits<Real>::infinity();
        Real best_left = 0.0;
        int best_segment = -1;
        int best_interval = -1;
        std::vector<std::pair<int, int>> best_reservations;
        const int maximum_radius = std::max(row_search_limit,
                                             static_cast<int>(db.rows.size()));
        for (int radius = 0; radius <= maximum_radius; ++radius) {
            bool tested = false;
            for (int direction : {-1, 1}) {
                if (radius == 0 && direction == 1) continue;
                const int r = base_row + direction * radius;
                if (r < 0 || r >= static_cast<int>(db.rows.size())) continue;
                const Row& row = db.rows[r];
                tested = true;
                if (node.height > row.height + 1.0e-6) {
                    Real left = 0.0;
                    std::vector<std::pair<int, int>> reservations;
                    if (!fit_multiline_cell(
                            db, segments_by_row, segments, r, node,
                            desired[id].x - 0.5 * node.width,
                            left, reservations)) continue;
                    const Real cost = std::abs(
                        left + 0.5 * node.width - desired[id].x) +
                        2.0 * std::abs(row.y + 0.5 * node.height - desired[id].y);
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_left = left;
                        best_segment = reservations.front().first;
                        best_interval = reservations.front().second;
                        best_reservations = std::move(reservations);
                    }
                    continue;
                }
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
                        best_reservations.clear();
                    }
                }
            }
            if (best_segment >= 0 && (radius >= row_search_limit ||
                2.0 * radius * db.rows.front().height > best_cost)) break;
            if (!tested && radius >= static_cast<int>(db.rows.size())) break;
        }
        if (best_segment < 0) {
            return false;
        }
        Segment& segment = segments[best_segment];
        const Row& row = db.rows[segment.row];
        node.x = best_left + 0.5 * node.width;
        node.y = row.y + 0.5 * node.height;
        if (best_reservations.empty()) {
            occupy(segment, best_interval, best_left, best_left + node.width);
            segment.cells.push_back(id);
        } else {
            // Occupy every covered row.  The tall cell itself is deliberately
            // omitted from the one-row detailed-placement segment, so later
            // reorder/Abacus passes keep this legal multi-row anchor fixed.
            for (auto [s, interval] : best_reservations)
                occupy(segments[s], interval, best_left, best_left + node.width);
        }
        ++placed;
        }
        std::cout << "[Greedy] placed=" << placed << " segments=" << segments.size() << '\n';
        return true;
    };

    std::vector<int> locality_order = db.movable_ids;
    std::stable_sort(locality_order.begin(), locality_order.end(), [&](int a, int b) {
        if (desired[a].x != desired[b].x) return desired[a].x < desired[b].x;
        if (db.nodes[a].width != db.nodes[b].width)
            return db.nodes[a].width > db.nodes[b].width;
        return desired[a].y < desired[b].y;
    });
    if (try_order(locality_order)) return;

    // Highly fragmented macro designs can exhaust the only intervals that
    // fit a wide cell even though total capacity is sufficient.  Reordering
    // by width first gives those cells first access to the larger intervals,
    // while retaining x/y as deterministic tie breakers.
    std::vector<int> width_order = db.movable_ids;
    std::stable_sort(width_order.begin(), width_order.end(), [&](int a, int b) {
        if (db.nodes[a].width != db.nodes[b].width)
            return db.nodes[a].width > db.nodes[b].width;
        if (desired[a].x != desired[b].x) return desired[a].x < desired[b].x;
        return desired[a].y < desired[b].y;
    });
    if (try_order(width_order)) {
        std::cout << "[GreedyFallback] width-first ordering accepted\n";
        return;
    }

    restore_movable_coordinates(db, original);
    throw std::runtime_error("greedy legalization cannot place all movable cells");
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

std::vector<int> minimum_cost_assignment(
    const std::vector<std::vector<Real>>& cost) {
    const int n = static_cast<int>(cost.size());
    std::vector<Real> u(n + 1, 0.0), v(n + 1, 0.0);
    std::vector<int> p(n + 1, 0), way(n + 1, 0);
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<Real> minv(n + 1, std::numeric_limits<Real>::infinity());
        std::vector<char> used(n + 1, false);
        do {
            used[j0] = true;
            const int i0 = p[j0];
            Real delta = std::numeric_limits<Real>::infinity();
            int j1 = 0;
            for (int j = 1; j <= n; ++j) {
                if (used[j]) continue;
                const Real current = cost[i0 - 1][j - 1] - u[i0] - v[j];
                if (current < minv[j]) {
                    minv[j] = current;
                    way[j] = j0;
                }
                if (minv[j] < delta) {
                    delta = minv[j];
                    j1 = j;
                }
            }
            for (int j = 0; j <= n; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            const int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0 != 0);
    }
    std::vector<int> assignment(n, 0);
    for (int j = 1; j <= n; ++j) assignment[p[j] - 1] = j - 1;
    return assignment;
}

std::size_t independent_set_matching(Database& db, int maximum_size,
                                     bool use_hungarian) {
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
            std::vector<std::vector<Real>> costs(
                cells.size(), std::vector<Real>(cells.size(), 0.0));
            for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
                Node& node = db.nodes[cells[i]];
                const Real old_x = node.x, old_y = node.y;
                for (int j = 0; j < static_cast<int>(slots.size()); ++j) {
                    node.x = slots[j].first;
                    node.y = slots[j].second;
                    costs[i][j] = affected_hpwl(db, incident[cells[i]]);
                }
                node.x = old_x;
                node.y = old_y;
            }
            std::vector<int> best_permutation;
            if (use_hungarian || cells.size() > 6) {
                best_permutation = minimum_cost_assignment(costs);
            } else {
                std::vector<int> permutation(cells.size());
                std::iota(permutation.begin(), permutation.end(), 0);
                Real best_cost = std::numeric_limits<Real>::infinity();
                best_permutation = permutation;
                do {
                    Real cost = 0.0;
                    for (int i = 0; i < static_cast<int>(cells.size()); ++i)
                        cost += costs[i][permutation[i]];
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_permutation = permutation;
                    }
                } while (std::next_permutation(permutation.begin(), permutation.end()));
            }
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

std::size_t cell_insertion(Database& db, std::vector<Segment>& segments,
                           int window, int passes) {
    const std::vector<std::vector<int>> incident = build_incident_nets(db);
    std::size_t accepted_total = 0;
    for (int pass = 0; pass < passes; ++pass) {
        std::size_t accepted = 0;
        for (Segment& segment : segments) {
            std::sort(segment.cells.begin(), segment.cells.end(), [&](int a, int b) {
                return db.nodes[a].x < db.nodes[b].x;
            });
            for (int source = 0; source < static_cast<int>(segment.cells.size()); ++source) {
                const int id = segment.cells[source];
                const Real target = net_target(db, id, incident[id]).first;
                auto position = std::lower_bound(
                    segment.cells.begin(), segment.cells.end(), target,
                    [&](int candidate, Real x) { return db.nodes[candidate].x < x; });
                int destination = static_cast<int>(position - segment.cells.begin());
                destination = std::clamp(destination, 0,
                    static_cast<int>(segment.cells.size()) - 1);
                destination = std::clamp(destination, source - window, source + window);
                if (destination == source) continue;
                const int begin = std::min(source, destination);
                const int end = std::max(source, destination) + 1;
                std::vector<int> original(segment.cells.begin() + begin,
                                          segment.cells.begin() + end);
                std::vector<int> candidate = original;
                const int local_source = source - begin;
                const int local_destination = destination - begin;
                const int moved = candidate[local_source];
                candidate.erase(candidate.begin() + local_source);
                candidate.insert(candidate.begin() + local_destination, moved);
                const Real left = db.nodes[original.front()].x -
                                  0.5 * db.nodes[original.front()].width;
                std::vector<Real> gaps(original.size() - 1, 0.0);
                for (int i = 0; i + 1 < static_cast<int>(original.size()); ++i) {
                    gaps[i] = (db.nodes[original[i + 1]].x -
                               0.5 * db.nodes[original[i + 1]].width) -
                              (db.nodes[original[i]].x +
                               0.5 * db.nodes[original[i]].width);
                }
                const std::vector<int> nets = union_incident(incident, original);
                const Real before = affected_hpwl(db, nets);
                assign_window(db, candidate, left, gaps);
                const Real after = affected_hpwl(db, nets);
                if (after + 1.0e-9 < before) {
                    std::copy(candidate.begin(), candidate.end(),
                              segment.cells.begin() + begin);
                    ++accepted;
                } else {
                    assign_window(db, original, left, gaps);
                }
            }
        }
        accepted_total += accepted;
        std::cout << "[CellInsertion] pass=" << pass
                  << " accepted=" << accepted << '\n';
        if (accepted == 0) break;
    }
    return accepted_total;
}

std::size_t projected_hpwl_refine(Database& db, std::vector<Segment>& segments,
                                  int passes, Real initial_step_sites,
                                  Real active_set_radius,
                                  Real active_set_power) {
    std::size_t accepted = 0;
    Real step_sites = initial_step_sites;
    for (int pass = 0; pass < passes; ++pass) {
        std::vector<Real> gx, gy;
        if (active_set_radius > 0.0) {
            exact_hpwl_active_set_direction(
                db, 100, active_set_radius, active_set_power, &gx, &gy);
        } else {
            exact_hpwl_subgradient(db, 100, &gx, &gy);
        }
        std::vector<Desired> desired(db.nodes.size());
        std::vector<std::pair<Real, Real>> original(db.nodes.size());
        for (int id : db.movable_ids) {
            original[id] = {db.nodes[id].x, db.nodes[id].y};
            const int row_index = nearest_row(
                db, db.nodes[id].y - 0.5 * db.nodes[id].height);
            const Real site = db.rows[row_index].site_spacing;
            const Real scale = step_sites * site /
                std::max<Real>(1.0, db.node_pin_weight[id]);
            desired[id] = {db.nodes[id].x - scale * gx[id], db.nodes[id].y};
        }
        const Real before = exact_hpwl(db);
        for (Segment& segment : segments) abacus_segment(db, segment, desired);
        const Real after = exact_hpwl(db);
        const LegalityResult legal = check_legality(db);
        if (legal.legal && after + 1.0e-9 < before) {
            ++accepted;
            std::cout << "[ProjectedHPWL] pass=" << pass
                      << " before=" << before << " after=" << after
                      << " step_sites=" << step_sites << '\n';
        } else {
            for (int id : db.movable_ids) {
                db.nodes[id].x = original[id].first;
                db.nodes[id].y = original[id].second;
            }
            step_sites *= 0.5;
            std::cout << "[ProjectedHPWL] pass=" << pass
                      << " rejected step_sites=" << step_sites << '\n';
            if (step_sites < 0.25) break;
        }
    }
    return accepted;
}

struct LegalBundleCut {
    Real hpwl = 0.0;
    std::vector<Real> x;
    std::vector<Real> gx;
};

void append_legal_bundle_cut(const Database& db, Real hpwl,
                             const std::vector<Real>& gx, int maximum_cuts,
                             std::vector<LegalBundleCut>& cuts) {
    LegalBundleCut cut;
    cut.hpwl = hpwl;
    cut.x.resize(db.movable_ids.size());
    cut.gx.resize(db.movable_ids.size());
    for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
        const int id = db.movable_ids[i];
        cut.x[i] = db.nodes[id].x;
        cut.gx[i] = gx[id];
    }
    if (static_cast<int>(cuts.size()) >= maximum_cuts) cuts.erase(cuts.begin());
    cuts.push_back(std::move(cut));
}

std::vector<Real> legal_bundle_direction(const Database& db,
                                         const std::vector<LegalBundleCut>& cuts,
                                         Real step_sites) {
    const int count = static_cast<int>(cuts.size());
    std::vector<Real> beta(count, 0.0), gram(count * count, 0.0);
    std::vector<Real> metric(db.movable_ids.size(), 0.0);
    for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
        const int id = db.movable_ids[i];
        const int row = nearest_row(db, db.nodes[id].y - 0.5 * db.nodes[id].height);
        metric[i] = db.rows[row].site_spacing /
                    std::max(1, db.node_pin_weight[id]);
    }
    for (int j = 0; j < count; ++j) {
        beta[j] = cuts[j].hpwl;
        for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
            const int id = db.movable_ids[i];
            beta[j] += cuts[j].gx[i] * (db.nodes[id].x - cuts[j].x[i]);
        }
        for (int k = 0; k <= j; ++k) {
            Real value = 0.0;
            for (std::size_t i = 0; i < db.movable_ids.size(); ++i)
                value += metric[i] * cuts[j].gx[i] * cuts[k].gx[i];
            gram[j * count + k] = gram[k * count + j] = value;
        }
    }
    std::vector<Real> alpha(count, 0.0);
    alpha.back() = 1.0;
    for (int iteration = 0; iteration < 30; ++iteration) {
        std::vector<Real> dual_gradient(count, 0.0);
        int best = 0;
        for (int j = 0; j < count; ++j) {
            Real product = 0.0;
            for (int k = 0; k < count; ++k)
                product += gram[j * count + k] * alpha[k];
            dual_gradient[j] = beta[j] - step_sites * product;
            if (dual_gradient[j] > dual_gradient[best]) best = j;
        }
        std::vector<Real> direction(count, 0.0);
        for (int j = 0; j < count; ++j) direction[j] = -alpha[j];
        direction[best] += 1.0;
        Real numerator = 0.0, curvature = 0.0;
        for (int j = 0; j < count; ++j) {
            numerator += dual_gradient[j] * direction[j];
            for (int k = 0; k < count; ++k)
                curvature += direction[j] * gram[j * count + k] * direction[k];
        }
        if (numerator <= 1.0e-7) break;
        const Real fraction = curvature > 1.0e-20
            ? std::min<Real>(1.0, numerator / (step_sites * curvature)) : 1.0;
        for (int j = 0; j < count; ++j) alpha[j] += fraction * direction[j];
    }
    std::vector<Real> direction(db.movable_ids.size(), 0.0);
    for (int j = 0; j < count; ++j) {
        for (std::size_t i = 0; i < direction.size(); ++i)
            direction[i] -= step_sites * metric[i] * alpha[j] * cuts[j].gx[i];
    }
    return direction;
}

std::size_t constrained_legal_bundle_refine(
    Database& db, std::vector<Segment>& segments, int passes,
    int maximum_cuts, Real initial_step_sites) {
    std::vector<LegalBundleCut> cuts;
    Real step_sites = initial_step_sites;
    std::size_t accepted = 0;
    for (int pass = 0; pass < passes; ++pass) {
        std::vector<Real> gx, gy;
        const Real before = exact_hpwl_subgradient(db, 100, &gx, &gy);
        append_legal_bundle_cut(db, before, gx, maximum_cuts, cuts);
        const std::vector<Real> direction =
            legal_bundle_direction(db, cuts, step_sites);
        const MovableCoordinates original = save_movable_coordinates(db);
        std::vector<Desired> desired(db.nodes.size());
        for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
            const int id = db.movable_ids[i];
            desired[id] = {db.nodes[id].x + direction[i], db.nodes[id].y};
        }
        for (Segment& segment : segments) abacus_segment(db, segment, desired);
        const Real after = exact_hpwl(db);
        const LegalityResult legal = check_legality(db);
        if (legal.legal && after + 1.0e-9 < before) {
            ++accepted;
            step_sites = std::min(initial_step_sites * 4.0, step_sites * 1.15);
            std::cout << "[LegalBundleSerious] pass=" << pass
                      << " before=" << before << " after=" << after
                      << " cuts=" << cuts.size()
                      << " step_sites=" << step_sites << '\n';
        } else {
            std::vector<Real> trial_gx, trial_gy;
            exact_hpwl_subgradient(db, 100, &trial_gx, &trial_gy);
            append_legal_bundle_cut(db, after, trial_gx, maximum_cuts, cuts);
            restore_movable_coordinates(db, original);
            step_sites *= 0.5;
            std::cout << "[LegalBundleNull] pass=" << pass
                      << " before=" << before << " trial=" << after
                      << " legal=" << (legal.legal ? 1 : 0)
                      << " cuts=" << cuts.size()
                      << " step_sites=" << step_sites << '\n';
            if (step_sites < 0.125) break;
        }
    }
    return accepted;
}

bool row_relegalization(Database& db, std::vector<Segment>& segments,
                        int passes, int row_search_limit) {
    bool improved = false;
    for (int pass = 0; pass < passes; ++pass) {
        const std::vector<std::vector<int>> incident = build_incident_nets(db);
        std::vector<Desired> desired(db.nodes.size());
        for (int id : db.movable_ids) {
            const auto target = net_target(db, id, incident[id]);
            const int row_index = nearest_row(
                db, db.nodes[id].y - 0.5 * db.nodes[id].height);
            const Real row_height = db.rows[row_index].height;
            desired[id].x = db.nodes[id].x;
            desired[id].y = db.nodes[id].y +
                std::clamp(target.second - db.nodes[id].y,
                           -row_height, row_height);
        }
        Database candidate = db;
        std::vector<Segment> candidate_segments = build_segments(candidate);
        try {
            greedy_legalize(candidate, desired, candidate_segments, row_search_limit);
            for (Segment& segment : candidate_segments)
                abacus_segment(candidate, segment, desired);
        } catch (const std::exception& error) {
            std::cout << "[RowRelegalize] pass=" << pass
                      << " rejected=" << error.what() << '\n';
            break;
        }
        const Real before = exact_hpwl(db);
        const Real after = exact_hpwl(candidate);
        const LegalityResult legal = check_legality(candidate);
        if (legal.legal && after + 1.0e-9 < before) {
            db = std::move(candidate);
            segments = std::move(candidate_segments);
            improved = true;
            std::cout << "[RowRelegalize] pass=" << pass
                      << " before=" << before << " after=" << after << '\n';
        } else {
            std::cout << "[RowRelegalize] pass=" << pass
                      << " no_improvement before=" << before
                      << " candidate=" << after << '\n';
            break;
        }
    }
    return improved;
}

bool overlaps(Real alo, Real ahi, Real blo, Real bhi, Real tolerance) {
    return std::min(ahi, bhi) - std::max(alo, blo) > tolerance;
}

}  // namespace

LegalizationProxy evaluate_legalization_proxy(const Database& db) {
    LegalizationProxy proxy;
    if (db.movable_ids.empty() || db.rows.empty()) return proxy;
    const std::vector<Segment> segments = build_segments(db);
    std::vector<std::vector<int>> segments_by_row(db.rows.size());
    for (int s = 0; s < static_cast<int>(segments.size()); ++s)
        segments_by_row[segments[s].row].push_back(s);
    std::vector<Real> demand(segments.size(), 0.0);
    Real total_width = 0.0;
    Real unassigned_width = 0.0;
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real bottom = node.y - 0.5 * node.height;
        const int row_index = nearest_row(db, bottom);
        const Row& row = db.rows[row_index];
        proxy.mean_row_distance += std::abs(bottom - row.y) /
            std::max<Real>(row.height, 1.0e-9);
        total_width += node.width;
        int best_segment = -1;
        Real best_distance = std::numeric_limits<Real>::infinity();
        const Real left = node.x - 0.5 * node.width;
        const Real right = node.x + 0.5 * node.width;
        for (int s : segments_by_row[row_index]) {
            const Real distance = left >= segments[s].lo && right <= segments[s].hi
                ? 0.0 : std::min(std::abs(left - segments[s].hi),
                                 std::abs(right - segments[s].lo));
            if (distance < best_distance) {
                best_distance = distance;
                best_segment = s;
            }
        }
        if (best_segment >= 0) demand[best_segment] += node.width;
        else unassigned_width += node.width;
    }
    proxy.mean_row_distance /= static_cast<Real>(db.movable_ids.size());
    Real excess = unassigned_width;
    for (int s = 0; s < static_cast<int>(segments.size()); ++s)
        excess += std::max<Real>(0.0, demand[s] - (segments[s].hi - segments[s].lo));
    proxy.segment_overflow = excess / std::max<Real>(total_width, 1.0);
    return proxy;
}

void compute_legalization_force(const Database& db, Real congestion_gain,
                                std::vector<Real>& grad_x,
                                std::vector<Real>& grad_y) {
    grad_x.assign(db.nodes.size(), 0.0);
    grad_y.assign(db.nodes.size(), 0.0);
    if (db.rows.empty() || db.movable_ids.empty()) return;

    const std::vector<Segment> segments = build_segments(db);
    std::vector<std::vector<int>> by_row(db.rows.size());
    for (int s = 0; s < static_cast<int>(segments.size()); ++s) {
        by_row[segments[s].row].push_back(s);
    }
    std::vector<Real> demand(segments.size(), 0.0);

    auto closest_segment = [&](const Node& node, int row_index,
                               Real& target_x, Real& distance) {
        const std::vector<int>& row_segments = by_row[row_index];
        int best = -1;
        distance = std::numeric_limits<Real>::infinity();
        if (row_segments.empty()) return best;
        auto position = std::lower_bound(
            row_segments.begin(), row_segments.end(), node.x,
            [&](int segment, Real x) { return segments[segment].hi < x; });
        const int center = static_cast<int>(position - row_segments.begin());
        for (int offset = -2; offset <= 2; ++offset) {
            const int index = center + offset;
            if (index < 0 || index >= static_cast<int>(row_segments.size())) continue;
            const int s = row_segments[index];
            const Real lo = segments[s].lo + 0.5 * node.width;
            const Real hi = segments[s].hi - 0.5 * node.width;
            if (hi < lo) continue;
            const Real candidate_x = std::clamp(node.x, lo, hi);
            const Real candidate_distance = std::abs(node.x - candidate_x);
            if (candidate_distance < distance) {
                best = s;
                distance = candidate_distance;
                target_x = candidate_x;
            }
        }
        return best;
    };

    // Estimate segment load from the current continuous placement.
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const int row = nearest_row(db, node.y - 0.5 * node.height);
        Real target_x = node.x;
        Real distance = 0.0;
        const int segment = closest_segment(node, row, target_x, distance);
        if (segment >= 0) demand[segment] += node.width;
    }

    // Compare only the nearest row and its two neighbors. This remains a
    // continuous attraction field; it does not assign a row or segment.
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real bottom = node.y - 0.5 * node.height;
        const int base_row = nearest_row(db, bottom);
        Real best_cost = std::numeric_limits<Real>::infinity();
        Real best_x = node.x;
        Real best_y = bottom;
        for (int row = std::max(0, base_row - 1);
             row <= std::min(static_cast<int>(db.rows.size()) - 1, base_row + 1);
             ++row) {
            if (node.height > db.rows[row].height + 1.0e-6) continue;
            Real target_x = node.x;
            Real x_distance = 0.0;
            const int segment = closest_segment(node, row, target_x, x_distance);
            if (segment < 0) continue;
            const Real capacity = std::max<Real>(
                segments[segment].hi - segments[segment].lo, 1.0);
            const Real congestion = std::max<Real>(
                0.0, demand[segment] / capacity - 1.0);
            const Real row_height = std::max<Real>(db.rows[row].height, 1.0e-9);
            const Real cost = std::abs(bottom - db.rows[row].y) / row_height +
                x_distance / row_height + congestion_gain * congestion;
            if (cost < best_cost) {
                best_cost = cost;
                best_x = target_x;
                best_y = db.rows[row].y;
            }
        }
        const Real row_height = std::max<Real>(
            db.rows[base_row].height, 1.0e-9);
        const Real pin_scale = std::max(1, db.node_pin_weight[id]);
        grad_x[id] = pin_scale * std::clamp(
            (node.x - best_x) / row_height, -1.0, 1.0);
        grad_y[id] = pin_scale * std::clamp(
            (bottom - best_y) / row_height, -1.0, 1.0);
    }
}

LegalityResult check_legality(const Database& db, Real tolerance) {
    LegalityResult result;
    std::vector<std::vector<int>> by_row(db.rows.size());
    const std::vector<std::vector<int>> obstacles_by_row =
        index_fixed_obstacles_by_row(db);
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real left = node.x - 0.5 * node.width;
        const Real right = node.x + 0.5 * node.width;
        const Real bottom = node.y - 0.5 * node.height;
        int row_index = nearest_row(db, bottom);
        const Row& row = db.rows[row_index];
        const int row_count = covered_row_count(db, row_index, node);
        const Real top = bottom + node.height;
        bool boundary_ok = row_count > 0 &&
                           left >= row.origin - tolerance &&
                           right <= row.xh() + tolerance &&
                           std::abs(bottom - row.y) <= tolerance &&
                           top <= db.rows[row_index + row_count - 1].y +
                                      db.rows[row_index + row_count - 1].height + tolerance;
        if (!boundary_ok) {
            ++result.boundary_errors;
            if (result.first_error.empty()) result.first_error = "boundary: " + node.name;
        }
        const Real site = (left - row.origin) / row.site_spacing;
        if (std::abs(site - std::round(site)) > tolerance) {
            ++result.alignment_errors;
            if (result.first_error.empty()) result.first_error = "site alignment: " + node.name;
        }
        if (row_count > 0) {
            for (int covered = row_index; covered < row_index + row_count; ++covered) {
                by_row[covered].push_back(id);
                for (int fixed_id : obstacles_by_row[covered]) {
                    const Node& fixed = db.nodes[fixed_id];
                    if (overlaps(left, right, fixed.x - 0.5 * fixed.width,
                                 fixed.x + 0.5 * fixed.width, tolerance) &&
                        overlaps(bottom, top, fixed.y - 0.5 * fixed.height,
                                 fixed.y + 0.5 * fixed.height, tolerance)) {
                        ++result.overlap_errors;
                        if (result.first_error.empty()) result.first_error =
                            "fixed overlap: " + node.name + " / " + fixed.name;
                        break;
                    }
                }
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

    if (config.row_relegalization_passes > 0) {
        row_relegalization(db, segments, config.row_relegalization_passes,
                           config.row_search_limit);
    }
    const int outer_rounds = std::max(1, config.detailed_outer_rounds);
    Real previous_round = exact_hpwl(db);
    for (int round = 0; round < outer_rounds; ++round) {
        // Cross-row swaps in the previous round invalidate the original
        // segment membership.  Rebuild it from the accepted legal placement
        // before applying any order-based operator again.
        segments = rebuild_segment_cells(db);
        if (config.run_detailed) {
            const MovableCoordinates before = save_movable_coordinates(db);
            detailed_reorder(db, segments, config.detailed_passes);
            legal = check_legality(db);
            if (!legal.legal) {
                restore_movable_coordinates(db, before);
                segments = rebuild_segment_cells(db);
                std::cout << "[DetailedGuard] adjacent rollback: "
                          << legal.first_error << '\n';
            }
        }
        save_legal_snapshot(db, config, "adjacent");
        if (config.run_dreamplace_detailed) {
            const MovableCoordinates before = save_movable_coordinates(db);
            k_reorder(db, segments, config.k_reorder_size,
                      config.detailed_passes);
            legal = check_legality(db);
            if (!legal.legal) {
                restore_movable_coordinates(db, before);
                segments = rebuild_segment_cells(db);
                std::cout << "[DetailedGuard] k-reorder rollback: "
                          << legal.first_error << '\n';
            }
        }
        if (config.cell_insertion_passes > 0) {
            const MovableCoordinates before = save_movable_coordinates(db);
            cell_insertion(db, segments, config.cell_insertion_window,
                           config.cell_insertion_passes);
            legal = check_legality(db);
            if (!legal.legal) {
                restore_movable_coordinates(db, before);
                segments = rebuild_segment_cells(db);
                std::cout << "[DetailedGuard] insertion rollback: "
                          << legal.first_error << '\n';
            }
        }
        result.hpwl_after_k_reorder = exact_hpwl(db);
        save_legal_snapshot(db, config, "k_reorder");
        legal = check_legality(db);
        if (!legal.legal) throw std::runtime_error(
            "K-Reorder/insertion produced illegal placement: " + legal.first_error);

        if (config.run_dreamplace_detailed) {
            const MovableCoordinates before = save_movable_coordinates(db);
            global_swap(db, config.global_swap_passes);
            legal = check_legality(db);
            if (!legal.legal) {
                restore_movable_coordinates(db, before);
                std::cout << "[DetailedGuard] global-swap rollback: "
                          << legal.first_error << '\n';
            }
        }
        result.hpwl_after_global_swap = exact_hpwl(db);
        save_legal_snapshot(db, config, "global_swap");
        legal = check_legality(db);
        if (!legal.legal) throw std::runtime_error(
            "global swap produced illegal placement: " + legal.first_error);

        if (config.run_dreamplace_detailed) {
            const MovableCoordinates before = save_movable_coordinates(db);
            independent_set_matching(db, config.independent_set_size,
                                     config.use_hungarian_matching);
            legal = check_legality(db);
            if (!legal.legal) {
                restore_movable_coordinates(db, before);
                std::cout << "[DetailedGuard] independent-set rollback: "
                          << legal.first_error << '\n';
            }
        }
        segments = rebuild_segment_cells(db);
        if (config.projected_subgradient_passes > 0) {
            projected_hpwl_refine(db, segments,
                                  config.projected_subgradient_passes,
                                  config.projected_step_sites,
                                  config.projected_active_set_radius,
                                  config.projected_active_set_power);
        }
        if (config.constrained_bundle_passes > 0) {
            constrained_legal_bundle_refine(
                db, segments, config.constrained_bundle_passes,
                config.constrained_bundle_size,
                config.constrained_bundle_step_sites);
        }
        result.hpwl_after_constrained_bundle = exact_hpwl(db);
        result.hpwl_after_independent_set = exact_hpwl(db);
        save_legal_snapshot(db, config, "independent_set");
        legal = check_legality(db);
        if (!legal.legal) throw std::runtime_error(
            "detailed round produced illegal placement: " + legal.first_error);
        const Real current_round = result.hpwl_after_independent_set;
        const Real relative = (previous_round - current_round) /
            std::max<Real>(previous_round, 1.0);
        std::cout << "[DetailedOuter] round=" << round
                  << " hpwl=" << current_round
                  << " relative_improvement=" << relative << '\n';
        if (round + 1 < outer_rounds &&
            relative < config.minimum_relative_improvement) break;
        previous_round = current_round;
    }
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
