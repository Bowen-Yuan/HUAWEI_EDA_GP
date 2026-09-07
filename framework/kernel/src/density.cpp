#include "epsilon_active/density.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ea {
namespace {

struct IntervalInfo {
    Real overlap = 0.0;
    Real derivative = 0.0;
};

Real active_weight(Real value, Real power) {
    if (value <= 0.0) return 0.0;
    if (power == 1.0) return value;
    if (power == 2.0) return value * value;
    if (power == 4.0) {
        const Real square = value * value;
        return square * square;
    }
    return std::pow(value, power);
}

Real active_min_variable_derivative(Real variable, Real constant,
                                    Real epsilon, Real power) {
    if (epsilon <= 0.0) {
        if (variable < constant) return 1.0;
        if (variable > constant) return 0.0;
        return 0.5;
    }
    const Real variable_score = variable <= constant
        ? 1.0
        : active_weight(1.0 - (variable - constant) / epsilon, power);
    const Real constant_score = constant <= variable
        ? 1.0
        : active_weight(1.0 - (constant - variable) / epsilon, power);
    return variable_score / std::max<Real>(variable_score + constant_score, 1.0e-30);
}

Real active_max_variable_derivative(Real variable, Real constant,
                                    Real epsilon, Real power) {
    if (epsilon <= 0.0) {
        if (variable > constant) return 1.0;
        if (variable < constant) return 0.0;
        return 0.5;
    }
    const Real variable_score = variable >= constant
        ? 1.0
        : active_weight(1.0 - (constant - variable) / epsilon, power);
    const Real constant_score = constant >= variable
        ? 1.0
        : active_weight(1.0 - (variable - constant) / epsilon, power);
    return variable_score / std::max<Real>(variable_score + constant_score, 1.0e-30);
}

IntervalInfo interval_overlap(Real rectangle_low, Real rectangle_high,
                              Real bin_low, Real bin_high,
                              Real epsilon, Real power) {
    IntervalInfo result;
    result.overlap = std::max<Real>(
        0.0, std::min(rectangle_high, bin_high) -
             std::max(rectangle_low, bin_low));
    const Real high_derivative = active_min_variable_derivative(
        rectangle_high, bin_high, epsilon, power);
    const Real low_derivative = active_max_variable_derivative(
        rectangle_low, bin_low, epsilon, power);
    result.derivative = high_derivative - low_derivative;
    if (result.overlap <= 0.0) {
        const bool epsilon_adjacent =
            rectangle_high >= bin_low - epsilon &&
            rectangle_low <= bin_high + epsilon;
        if (!epsilon_adjacent) result.derivative = 0.0;
    }
    return result;
}

int lower_bin(Real coordinate, Real origin, Real pitch, int count) {
    return std::clamp(static_cast<int>(std::floor((coordinate - origin) / pitch)),
                      0, count - 1);
}

}  // namespace

std::vector<Real> prolongate_bin_field(
    const std::vector<Real>& source, int source_bins_x, int source_bins_y,
    int target_bins_x, int target_bins_y) {
    if (source_bins_x <= 0 || source_bins_y <= 0 || target_bins_x <= 0 ||
        target_bins_y <= 0 ||
        source.size() != static_cast<std::size_t>(source_bins_x) *
                             source_bins_y) {
        throw std::invalid_argument("invalid bin field dimensions");
    }
    std::vector<Real> target(
        static_cast<std::size_t>(target_bins_x) * target_bins_y, 0.0);
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < target_bins_y; ++y) {
        const int source_y = std::min(
            source_bins_y - 1,
            static_cast<int>((static_cast<long long>(y) * source_bins_y) /
                             target_bins_y));
        for (int x = 0; x < target_bins_x; ++x) {
            const int source_x = std::min(
                source_bins_x - 1,
                static_cast<int>((static_cast<long long>(x) * source_bins_x) /
                                 target_bins_x));
            target[static_cast<std::size_t>(y) * target_bins_x + x] =
                source[static_cast<std::size_t>(source_y) * source_bins_x +
                       source_x];
        }
    }
    return target;
}

ExactOverlapDensity::ExactOverlapDensity(const Database& db, int bins_x,
                                         int bins_y, Real target_density)
    : db_(db), bins_x_(bins_x), bins_y_(bins_y),
      target_density_(target_density),
      bin_width_((db.xh - db.xl) / bins_x),
      bin_height_((db.yh - db.yl) / bins_y),
      bin_area_(bin_width_ * bin_height_),
      fixed_occupancy_(static_cast<std::size_t>(bins_x) * bins_y, 0.0),
      occupancy_(fixed_occupancy_.size(), 0.0),
      excess_density_(fixed_occupancy_.size(), 0.0) {
    if (bins_x <= 0 || bins_y <= 0 || target_density <= 0.0 ||
        target_density > 1.0 || bin_width_ <= 0.0 || bin_height_ <= 0.0) {
        throw std::invalid_argument("invalid exact-overlap density grid");
    }
    for (int id : db_.fixed_ids) {
        const Node& node = db_.nodes[id];
        if (!node.terminal_ni) add_rectangle(fixed_occupancy_, node);
    }
}

void ExactOverlapDensity::add_rectangle(std::vector<Real>& map,
                                        const Node& node) const {
    const Real left = node.x - 0.5 * node.width;
    const Real right = node.x + 0.5 * node.width;
    const Real bottom = node.y - 0.5 * node.height;
    const Real top = node.y + 0.5 * node.height;
    if (right <= db_.xl || left >= db_.xh || top <= db_.yl || bottom >= db_.yh) return;
    const int x0 = lower_bin(left, db_.xl, bin_width_, bins_x_);
    const int x1 = lower_bin(std::nextafter(right, left), db_.xl, bin_width_, bins_x_);
    const int y0 = lower_bin(bottom, db_.yl, bin_height_, bins_y_);
    const int y1 = lower_bin(std::nextafter(top, bottom), db_.yl, bin_height_, bins_y_);
    for (int by = y0; by <= y1; ++by) {
        const Real bin_bottom = db_.yl + by * bin_height_;
        const Real overlap_y = std::max<Real>(
            0.0, std::min(top, bin_bottom + bin_height_) - std::max(bottom, bin_bottom));
        for (int bx = x0; bx <= x1; ++bx) {
            const Real bin_left = db_.xl + bx * bin_width_;
            const Real overlap_x = std::max<Real>(
                0.0, std::min(right, bin_left + bin_width_) - std::max(left, bin_left));
            map[static_cast<std::size_t>(by) * bins_x_ + bx] += overlap_x * overlap_y;
        }
    }
}

void ExactOverlapDensity::append_rectangle_changes(
    const Node& node, Real x, Real y, Real sign,
    std::vector<std::pair<int, Real>>& changes) const {
    const Real left = x - 0.5 * node.width;
    const Real right = x + 0.5 * node.width;
    const Real bottom = y - 0.5 * node.height;
    const Real top = y + 0.5 * node.height;
    if (right <= db_.xl || left >= db_.xh || top <= db_.yl || bottom >= db_.yh) return;
    const int x0 = lower_bin(left, db_.xl, bin_width_, bins_x_);
    const int x1 = lower_bin(std::nextafter(right, left), db_.xl, bin_width_, bins_x_);
    const int y0 = lower_bin(bottom, db_.yl, bin_height_, bins_y_);
    const int y1 = lower_bin(std::nextafter(top, bottom), db_.yl, bin_height_, bins_y_);
    for (int by = y0; by <= y1; ++by) {
        const Real bin_bottom = db_.yl + by * bin_height_;
        const Real overlap_y = std::max<Real>(
            0.0, std::min(top, bin_bottom + bin_height_) - std::max(bottom, bin_bottom));
        for (int bx = x0; bx <= x1; ++bx) {
            const Real bin_left = db_.xl + bx * bin_width_;
            const Real overlap_x = std::max<Real>(
                0.0, std::min(right, bin_left + bin_width_) - std::max(left, bin_left));
            if (overlap_x > 0.0 && overlap_y > 0.0) {
                changes.emplace_back(by * bins_x_ + bx,
                                     sign * overlap_x * overlap_y);
            }
        }
    }
}

DensityMove ExactOverlapDensity::evaluate_move(int node_id, Real new_x,
                                                Real new_y) const {
    if (node_id < 0 || node_id >= static_cast<int>(db_.nodes.size()) ||
        db_.nodes[node_id].fixed) {
        throw std::invalid_argument("density move requires a movable node");
    }
    DensityMove result;
    const Node& node = db_.nodes[node_id];
    append_rectangle_changes(node, node.x, node.y, -1.0, result.area_changes);
    append_rectangle_changes(node, new_x, new_y, 1.0, result.area_changes);
    std::sort(result.area_changes.begin(), result.area_changes.end());
    std::size_t write = 0;
    for (std::size_t read = 0; read < result.area_changes.size();) {
        const int bin = result.area_changes[read].first;
        Real delta = 0.0;
        while (read < result.area_changes.size() &&
               result.area_changes[read].first == bin) {
            delta += result.area_changes[read].second;
            ++read;
        }
        if (std::abs(delta) > 1.0e-14) {
            result.area_changes[write++] = {bin, delta};
        }
    }
    result.area_changes.resize(write);
    const Real capacity = target_density_ * bin_area_;
    for (const auto& [bin, delta] : result.area_changes) {
        const Real old_excess = std::max<Real>(occupancy_[bin] - capacity, 0.0);
        const Real new_excess = std::max<Real>(occupancy_[bin] + delta - capacity, 0.0);
        result.overflow_area_delta += new_excess - old_excess;
        result.energy_delta += 0.5 / bin_area_ *
            (new_excess * new_excess - old_excess * old_excess);
    }
    return result;
}

DensityMove ExactOverlapDensity::evaluate_group_move(
    const std::vector<DensityNodeMove>& moves) const {
    DensityMove result;
    for (const DensityNodeMove& move : moves) {
        if (move.node_id < 0 ||
            move.node_id >= static_cast<int>(db_.nodes.size()) ||
            db_.nodes[move.node_id].fixed) {
            throw std::invalid_argument(
                "density group move requires movable nodes");
        }
        const Node& node = db_.nodes[move.node_id];
        append_rectangle_changes(
            node, node.x, node.y, -1.0, result.area_changes);
        append_rectangle_changes(
            node, move.x, move.y, 1.0, result.area_changes);
    }
    std::sort(result.area_changes.begin(), result.area_changes.end());
    std::size_t write = 0;
    for (std::size_t read = 0; read < result.area_changes.size();) {
        const int bin = result.area_changes[read].first;
        Real delta = 0.0;
        while (read < result.area_changes.size() &&
               result.area_changes[read].first == bin) {
            delta += result.area_changes[read].second;
            ++read;
        }
        if (std::abs(delta) > 1.0e-14) {
            result.area_changes[write++] = {bin, delta};
        }
    }
    result.area_changes.resize(write);
    const Real capacity = target_density_ * bin_area_;
    for (const auto& [bin, delta] : result.area_changes) {
        const Real old_excess = std::max<Real>(occupancy_[bin] - capacity, 0.0);
        const Real new_excess =
            std::max<Real>(occupancy_[bin] + delta - capacity, 0.0);
        result.overflow_area_delta += new_excess - old_excess;
        result.energy_delta += 0.5 / bin_area_ *
            (new_excess * new_excess - old_excess * old_excess);
    }
    return result;
}

void ExactOverlapDensity::commit_move(const DensityMove& move) {
    for (const auto& [bin, delta] : move.area_changes) occupancy_[bin] += delta;
}

DensityMetrics ExactOverlapDensity::evaluate(Real epsilon, Real active_power,
                                              std::vector<Real>* grad_x,
                                              std::vector<Real>* grad_y) {
    if (epsilon < 0.0 || active_power <= 0.0) {
        throw std::invalid_argument("invalid density active-set parameters");
    }
    occupancy_ = fixed_occupancy_;
    #pragma omp parallel for schedule(static)
    for (int movable = 0; movable < static_cast<int>(db_.movable_ids.size()); ++movable) {
        const Node& node = db_.nodes[db_.movable_ids[movable]];
        const Real left = node.x - 0.5 * node.width;
        const Real right = node.x + 0.5 * node.width;
        const Real bottom = node.y - 0.5 * node.height;
        const Real top = node.y + 0.5 * node.height;
        const int x0 = lower_bin(left, db_.xl, bin_width_, bins_x_);
        const int x1 = lower_bin(std::nextafter(right, left), db_.xl, bin_width_, bins_x_);
        const int y0 = lower_bin(bottom, db_.yl, bin_height_, bins_y_);
        const int y1 = lower_bin(std::nextafter(top, bottom), db_.yl, bin_height_, bins_y_);
        for (int by = y0; by <= y1; ++by) {
            const Real bin_bottom = db_.yl + by * bin_height_;
            const Real overlap_y = std::max<Real>(
                0.0, std::min(top, bin_bottom + bin_height_) - std::max(bottom, bin_bottom));
            for (int bx = x0; bx <= x1; ++bx) {
                const Real bin_left = db_.xl + bx * bin_width_;
                const Real overlap_x = std::max<Real>(
                    0.0, std::min(right, bin_left + bin_width_) - std::max(left, bin_left));
                const Real area = overlap_x * overlap_y;
                const std::size_t index = static_cast<std::size_t>(by) * bins_x_ + bx;
                #pragma omp atomic update
                occupancy_[index] += area;
            }
        }
    }

    DensityMetrics metrics;
    Real excess_area = 0.0;
    Real maximum_density = 0.0;
    Real energy = 0.0;
    #pragma omp parallel for reduction(+:excess_area,energy) reduction(max:maximum_density) schedule(static)
    for (int index = 0; index < static_cast<int>(occupancy_.size()); ++index) {
        const Real density = occupancy_[index] / bin_area_;
        const Real excess = std::max<Real>(density - target_density_, 0.0);
        excess_density_[index] = excess;
        excess_area += excess * bin_area_;
        energy += 0.5 * bin_area_ * excess * excess;
        maximum_density = std::max(maximum_density, density);
    }
    metrics.energy = energy;
    metrics.max_density = maximum_density;
    metrics.overflow = excess_area / std::max<Real>(db_.movable_area, 1.0e-30);

    if (grad_x) grad_x->assign(db_.nodes.size(), 0.0);
    if (grad_y) grad_y->assign(db_.nodes.size(), 0.0);
    if (!grad_x && !grad_y) return metrics;

    #pragma omp parallel for schedule(static)
    for (int movable = 0; movable < static_cast<int>(db_.movable_ids.size()); ++movable) {
        const int node_id = db_.movable_ids[movable];
        const Node& node = db_.nodes[node_id];
        const Real left = node.x - 0.5 * node.width;
        const Real right = node.x + 0.5 * node.width;
        const Real bottom = node.y - 0.5 * node.height;
        const Real top = node.y + 0.5 * node.height;
        const int x0 = lower_bin(left - epsilon, db_.xl, bin_width_, bins_x_);
        const int x1 = lower_bin(std::nextafter(right + epsilon, left),
                                 db_.xl, bin_width_, bins_x_);
        const int y0 = lower_bin(bottom - epsilon, db_.yl, bin_height_, bins_y_);
        const int y1 = lower_bin(std::nextafter(top + epsilon, bottom),
                                 db_.yl, bin_height_, bins_y_);
        Real gx = 0.0;
        Real gy = 0.0;
        for (int by = y0; by <= y1; ++by) {
            const Real bin_bottom = db_.yl + by * bin_height_;
            const IntervalInfo iy = interval_overlap(
                bottom, top, bin_bottom, bin_bottom + bin_height_, epsilon, active_power);
            for (int bx = x0; bx <= x1; ++bx) {
                const std::size_t index = static_cast<std::size_t>(by) * bins_x_ + bx;
                const Real coefficient = excess_density_[index];
                if (coefficient <= 0.0) continue;
                const Real bin_left = db_.xl + bx * bin_width_;
                const IntervalInfo ix = interval_overlap(
                    left, right, bin_left, bin_left + bin_width_, epsilon, active_power);
                gx += coefficient * ix.derivative * iy.overlap;
                gy += coefficient * iy.derivative * ix.overlap;
            }
        }
        if (grad_x) (*grad_x)[node_id] = gx;
        if (grad_y) (*grad_y)[node_id] = gy;
    }
    return metrics;
}

DensityMetrics ExactOverlapDensity::evaluate_with_prices(
    const std::vector<Real>& prices, Real epsilon, Real active_power,
    std::vector<Real>* grad_x, std::vector<Real>* grad_y) {
    if (prices.size() != occupancy_.size()) {
        throw std::invalid_argument("regional price field has wrong size");
    }
    // Rebuild occupancy and exact excess first. This deliberately delegates
    // the metric calculation to the unweighted oracle; prices never alter
    // the reported density or overflow.
    const DensityMetrics metrics = evaluate(epsilon, active_power, nullptr, nullptr);
    if (!grad_x && !grad_y) return metrics;
    if (grad_x) grad_x->assign(db_.nodes.size(), 0.0);
    if (grad_y) grad_y->assign(db_.nodes.size(), 0.0);
    #pragma omp parallel for schedule(static)
    for (int movable = 0; movable < static_cast<int>(db_.movable_ids.size()); ++movable) {
        const int node_id = db_.movable_ids[movable];
        const Node& node = db_.nodes[node_id];
        const Real left = node.x - 0.5 * node.width;
        const Real right = node.x + 0.5 * node.width;
        const Real bottom = node.y - 0.5 * node.height;
        const Real top = node.y + 0.5 * node.height;
        const int x0 = lower_bin(left - epsilon, db_.xl, bin_width_, bins_x_);
        const int x1 = lower_bin(std::nextafter(right + epsilon, left),
                                 db_.xl, bin_width_, bins_x_);
        const int y0 = lower_bin(bottom - epsilon, db_.yl, bin_height_, bins_y_);
        const int y1 = lower_bin(std::nextafter(top + epsilon, bottom),
                                 db_.yl, bin_height_, bins_y_);
        Real gx = 0.0;
        Real gy = 0.0;
        for (int by = y0; by <= y1; ++by) {
            const Real bin_bottom = db_.yl + by * bin_height_;
            const IntervalInfo iy = interval_overlap(
                bottom, top, bin_bottom, bin_bottom + bin_height_, epsilon, active_power);
            for (int bx = x0; bx <= x1; ++bx) {
                const int index = by * bins_x_ + bx;
                const Real coefficient = excess_density_[index] * prices[index];
                if (coefficient <= 0.0) continue;
                const Real bin_left = db_.xl + bx * bin_width_;
                const IntervalInfo ix = interval_overlap(
                    left, right, bin_left, bin_left + bin_width_, epsilon, active_power);
                gx += coefficient * ix.derivative * iy.overlap;
                gy += coefficient * iy.derivative * ix.overlap;
            }
        }
        if (grad_x) (*grad_x)[node_id] = gx;
        if (grad_y) (*grad_y)[node_id] = gy;
    }
    return metrics;
}

}  // namespace ea
