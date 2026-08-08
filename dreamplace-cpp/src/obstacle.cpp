#include "obstacle.h"

#include <algorithm>
#include <cmath>

namespace dpcpp {

FixedObstacleField::FixedObstacleField(const Database& db, int bins_x,
                                       int bins_y)
    : nx_(std::max(1, bins_x)), ny_(std::max(1, bins_y)),
      xl_(db.xl), yl_(db.yl), xh_(db.xh), yh_(db.yh) {
    bin_w_ = (xh_ - xl_) / nx_;
    bin_h_ = (yh_ - yl_) / ny_;
    candidates_.resize(nx_ * ny_);

    Real max_half_width = 0.0;
    Real max_half_height = 0.0;
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        movable_area_ += node.area();
        max_half_width = std::max(max_half_width, 0.5 * node.width);
        max_half_height = std::max(max_half_height, 0.5 * node.height);
    }
    for (int fixed_id : db.fixed_ids) {
        const Node& fixed = db.nodes[fixed_id];
        if (fixed.terminal_ni) continue;
        const Real left = fixed.x - 0.5 * fixed.width - max_half_width;
        const Real right = fixed.x + 0.5 * fixed.width + max_half_width;
        const Real bottom = fixed.y - 0.5 * fixed.height - max_half_height;
        const Real top = fixed.y + 0.5 * fixed.height + max_half_height;
        const int bx0 = std::clamp(
            static_cast<int>(std::floor((left - xl_) / bin_w_)), 0, nx_ - 1);
        const int bx1 = std::clamp(
            static_cast<int>(std::floor((right - xl_) / bin_w_)), 0, nx_ - 1);
        const int by0 = std::clamp(
            static_cast<int>(std::floor((bottom - yl_) / bin_h_)), 0, ny_ - 1);
        const int by1 = std::clamp(
            static_cast<int>(std::floor((top - yl_) / bin_h_)), 0, ny_ - 1);
        for (int by = by0; by <= by1; ++by) {
            for (int bx = bx0; bx <= bx1; ++bx) {
                candidates_[by * nx_ + bx].push_back(fixed_id);
            }
        }
    }
}

ObstacleResult FixedObstacleField::compute(const Database& db,
                                           std::vector<Real>* grad_x,
                                           std::vector<Real>* grad_y) const {
    if (grad_x) grad_x->assign(db.nodes.size(), 0.0);
    if (grad_y) grad_y->assign(db.nodes.size(), 0.0);
    ObstacleResult result;
    const Real mid_x = 0.5 * (xl_ + xh_);
    const Real mid_y = 0.5 * (yl_ + yh_);

    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const int bx = std::clamp(
            static_cast<int>(std::floor((node.x - xl_) / bin_w_)), 0, nx_ - 1);
        const int by = std::clamp(
            static_cast<int>(std::floor((node.y - yl_) / bin_h_)), 0, ny_ - 1);
        bool overlaps_any = false;
        for (int fixed_id : candidates_[by * nx_ + bx]) {
            const Node& fixed = db.nodes[fixed_id];
            const Real dx = node.x - fixed.x;
            const Real dy = node.y - fixed.y;
            const Real penetration_x = 0.5 * (node.width + fixed.width) - std::abs(dx);
            const Real penetration_y = 0.5 * (node.height + fixed.height) - std::abs(dy);
            if (penetration_x <= 0.0 || penetration_y <= 0.0) continue;

            overlaps_any = true;
            result.overlap_area += penetration_x * penetration_y;
            // Squared distance to the complement of the fixed rectangle
            // expanded by the movable half-size. The closest escape axis is
            // used, so gradient descent moves the cell out without snapping.
            if (penetration_x <= penetration_y) {
                const Real direction = std::abs(dx) > 1.0e-12
                    ? (dx > 0.0 ? 1.0 : -1.0)
                    : (fixed.x <= mid_x ? 1.0 : -1.0);
                result.energy += 0.5 * penetration_x * penetration_x;
                if (grad_x) (*grad_x)[id] -= direction * penetration_x;
            } else {
                const Real direction = std::abs(dy) > 1.0e-12
                    ? (dy > 0.0 ? 1.0 : -1.0)
                    : (fixed.y <= mid_y ? 1.0 : -1.0);
                result.energy += 0.5 * penetration_y * penetration_y;
                if (grad_y) (*grad_y)[id] -= direction * penetration_y;
            }
        }
        if (overlaps_any) ++result.overlapping_cells;
    }
    result.overlap_ratio = result.overlap_area /
        std::max<Real>(1.0, movable_area_);
    return result;
}

}  // namespace dpcpp
