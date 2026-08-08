#pragma once

#include "types.h"

#include <vector>

namespace dpcpp {

struct ObstacleResult {
    Real energy = 0.0;
    Real overlap_area = 0.0;
    Real overlap_ratio = 0.0;
    int overlapping_cells = 0;
};

// Continuous fixed-obstacle model used only by the opt-in progressive flow.
// It indexes fixed rectangles spatially and differentiates the squared
// distance that a movable center must travel to leave each expanded macro.
class FixedObstacleField {
public:
    explicit FixedObstacleField(const Database& db, int bins_x = 64,
                                int bins_y = 64);

    ObstacleResult compute(const Database& db, std::vector<Real>* grad_x,
                           std::vector<Real>* grad_y) const;

private:
    int nx_ = 0;
    int ny_ = 0;
    Real xl_ = 0.0;
    Real yl_ = 0.0;
    Real xh_ = 0.0;
    Real yh_ = 0.0;
    Real bin_w_ = 0.0;
    Real bin_h_ = 0.0;
    Real movable_area_ = 0.0;
    std::vector<std::vector<int>> candidates_;
};

}  // namespace dpcpp
