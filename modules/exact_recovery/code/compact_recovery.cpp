#include "epsilon_active/compact_recovery.hpp"

#include "epsilon_active/bookshelf.hpp"
#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ea {
namespace {

struct Point {
    Real x = 0.0;
    Real y = 0.0;
};

Point movable_centroid(const Database& db) {
    Real weight = 0.0;
    Point centroid;
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real area = std::max<Real>(node.width * node.height, 1.0);
        weight += area;
        centroid.x += area * node.x;
        centroid.y += area * node.y;
    }
    if (weight <= 0.0) {
        return {0.5 * (db.xl + db.xh), 0.5 * (db.yl + db.yh)};
    }
    centroid.x /= weight;
    centroid.y /= weight;
    return centroid;
}

Real compactness(const Database& db, const Point& centroid) {
    Real total = 0.0;
    Real weight = 0.0;
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real area = std::max<Real>(node.width * node.height, 1.0);
        const Real dx = node.x - centroid.x;
        const Real dy = node.y - centroid.y;
        total += area * (dx * dx + dy * dy);
        weight += area;
    }
    return total / std::max<Real>(weight, 1.0);
}

}  // namespace

CompactRecoveryStats compact_support_contraction(
    Database& db, ExactOverlapDensity& density,
    const CompactRecoveryConfig& config) {
    if (config.sweeps < 0 || config.line_search_steps <= 0 ||
        config.max_step_bins <= 0.0 || config.overflow_cap >= 1.0 ||
        config.hpwl_tolerance < 0.0 || config.snapshot_every < 0) {
        throw std::invalid_argument("invalid compact recovery configuration");
    }
    CompactRecoveryStats stats;
    if (config.sweeps == 0) return stats;
    const auto started = std::chrono::steady_clock::now();
    ExactHpwl hpwl(db);
    DensityMetrics current_density = density.evaluate(0.0, 1.0, nullptr, nullptr);
    Real current_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
    stats.objective_evaluations += 2;
    const Point center = movable_centroid(db);
    stats.initial_hpwl = current_hpwl;
    stats.initial_overflow = current_density.overflow;
    stats.initial_compactness = compactness(db, center);
    const Real cap = config.overflow_cap >= 0.0
        ? config.overflow_cap : current_density.overflow;
    const Real overflow_tolerance = 1.0e-10 * density.bin_area();
    std::vector<unsigned char> anchored(db.nodes.size(), 0);
    if (config.fixed_terminal_aware) {
        for (const Net& net : db.nets) {
            bool has_fixed = false;
            for (std::size_t p = net.pin_begin;
                 p < net.pin_begin + net.pin_count; ++p) {
                has_fixed = has_fixed || db.nodes[db.pins[p].node].fixed;
            }
            if (!has_fixed) continue;
            for (std::size_t p = net.pin_begin;
                 p < net.pin_begin + net.pin_count; ++p) {
                const int id = db.pins[p].node;
                if (!db.nodes[id].fixed) anchored[id] = 1;
            }
        }
    }

    for (int sweep = 0; sweep < config.sweeps; ++sweep) {
        const Point target = movable_centroid(db);
        const Real before_compactness = compactness(db, target);
        bool accepted = false;
        Real best_hpwl = current_hpwl;
        Real best_compactness = before_compactness;
        std::vector<DensityNodeMove> best_moves;

        const Real base = config.max_step_bins /
            static_cast<Real>(std::max(density.bins_x(), density.bins_y()));
        for (int line = 0; line < config.line_search_steps; ++line) {
            const Real scale = std::ldexp(base, -line);
            const std::vector<int> modes = config.axis_separated
                ? std::vector<int>{0, 1, 2} : std::vector<int>{0};
            for (int mode : modes) {
                ++stats.candidates;
                std::vector<DensityNodeMove> moves;
                moves.reserve(db.movable_ids.size());
                for (int id : db.movable_ids) {
                    const Node& node = db.nodes[id];
                    if (anchored[id]) {
                        moves.push_back({id, node.x, node.y});
                        continue;
                    }
                    const Real alpha = std::clamp(scale, 0.0, 0.5);
                    Real nx = node.x + alpha * (target.x - node.x);
                    Real ny = node.y + alpha * (target.y - node.y);
                    if (mode == 1) ny = node.y;
                    if (mode == 2) nx = node.x;
                    nx = std::clamp(nx, db.xl + 0.5 * node.width,
                                    db.xh - 0.5 * node.width);
                    ny = std::clamp(ny, db.yl + 0.5 * node.height,
                                    db.yh - 0.5 * node.height);
                    moves.push_back({id, nx, ny});
                }
                DensityMove density_move = density.evaluate_group_move(moves);
                const Real candidate_overflow = current_density.overflow +
                    density_move.overflow_area_delta /
                    std::max<Real>(db.movable_area, 1.0e-30);
                if (current_density.overflow > cap + overflow_tolerance) {
                    if (candidate_overflow > current_density.overflow +
                        overflow_tolerance / std::max<Real>(db.movable_area, 1.0)) {
                        continue;
                    }
                } else if (candidate_overflow > cap +
                           overflow_tolerance / std::max<Real>(db.movable_area, 1.0)) {
                    continue;
                }
                std::vector<Point> old_positions;
                old_positions.reserve(db.movable_ids.size());
                for (int id : db.movable_ids) {
                    old_positions.push_back({db.nodes[id].x, db.nodes[id].y});
                }
                for (const DensityNodeMove& move : moves) {
                    db.nodes[move.node_id].x = move.x;
                    db.nodes[move.node_id].y = move.y;
                }
                const Real candidate_hpwl = hpwl.evaluate(
                    0.0, 1.0, -1, nullptr, nullptr);
                const Real candidate_compactness = compactness(db, target);
                ++stats.objective_evaluations;
                for (std::size_t i = 0; i < db.movable_ids.size(); ++i) {
                    const int id = db.movable_ids[i];
                    db.nodes[id].x = old_positions[i].x;
                    db.nodes[id].y = old_positions[i].y;
                }
                current_density = density.evaluate(0.0, 1.0, nullptr, nullptr);
                current_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
                ++stats.objective_evaluations;
                if (candidate_hpwl + config.hpwl_tolerance < best_hpwl &&
                    (candidate_compactness < best_compactness ||
                     best_moves.empty())) {
                    accepted = true;
                    best_hpwl = candidate_hpwl;
                    best_compactness = candidate_compactness;
                    best_moves = std::move(moves);
                }
            }
        }
        if (accepted) {
            // Apply the selected contraction from the original state. The
            // current state was restored after every trial.
            DensityMove committed = density.evaluate_group_move(best_moves);
            density.commit_move(committed);
            for (const DensityNodeMove& move : best_moves) {
                db.nodes[move.node_id].x = move.x;
                db.nodes[move.node_id].y = move.y;
            }
            current_density = density.evaluate(0.0, 1.0, nullptr, nullptr);
            current_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
            ++stats.accepted_candidates;
            ++stats.moves;
        }
        stats.trajectory.push_back({sweep, accepted ? 1 : 0, current_hpwl,
                                    current_density.overflow,
                                    compactness(db, target)});
        ++stats.sweeps;
        if (config.snapshot_every > 0 &&
            (sweep + 1) % config.snapshot_every == 0 &&
            !config.output_dir.empty()) {
            std::filesystem::create_directories(config.output_dir / "snapshots");
            write_bookshelf_placement(db, config.output_dir / "snapshots" /
                ("compact_" + std::to_string(sweep + 1) + ".pl"));
        }
        if (!accepted) break;
    }
    stats.final_hpwl = current_hpwl;
    stats.final_overflow = current_density.overflow;
    stats.final_compactness = compactness(db, movable_centroid(db));
    stats.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    return stats;
}

}  // namespace ea
