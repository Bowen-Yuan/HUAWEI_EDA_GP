#pragma once

#include "epsilon_active/types.hpp"

#include <vector>
#include <utility>

namespace ea {

std::vector<Real> prolongate_bin_field(
    const std::vector<Real>& source, int source_bins_x, int source_bins_y,
    int target_bins_x, int target_bins_y);

struct DensityMove {
    Real overflow_area_delta = 0.0;
    Real energy_delta = 0.0;
    std::vector<std::pair<int, Real>> area_changes;
};

struct DensityNodeMove {
    int node_id = -1;
    Real x = 0.0;
    Real y = 0.0;
};

class ExactOverlapDensity {
public:
    ExactOverlapDensity(const Database& db, int bins_x, int bins_y,
                        Real target_density);

    DensityMetrics evaluate(Real epsilon, Real active_power,
                            std::vector<Real>* grad_x,
                            std::vector<Real>* grad_y);
    // Exact overlap evaluation with a bin-wise algorithmic price field. The
    // price only selects a search direction; occupancy and overflow remain
    // the unweighted exact rectangle/bin quantities.
    DensityMetrics evaluate_with_prices(
        const std::vector<Real>& prices, Real epsilon, Real active_power,
        std::vector<Real>* grad_x, std::vector<Real>* grad_y);
    DensityMove evaluate_move(int node_id, Real new_x, Real new_y) const;
    DensityMove evaluate_group_move(
        const std::vector<DensityNodeMove>& moves) const;
    void commit_move(const DensityMove& move);

    int bins_x() const noexcept { return bins_x_; }
    int bins_y() const noexcept { return bins_y_; }
    Real bin_width() const noexcept { return bin_width_; }
    Real bin_height() const noexcept { return bin_height_; }
    Real bin_area() const noexcept { return bin_area_; }
    Real target_density() const noexcept { return target_density_; }
    const std::vector<Real>& fixed_occupancy() const noexcept {
        return fixed_occupancy_;
    }
    const std::vector<Real>& occupancy() const noexcept { return occupancy_; }

private:
    const Database& db_;
    int bins_x_;
    int bins_y_;
    Real target_density_;
    Real bin_width_;
    Real bin_height_;
    Real bin_area_;
    std::vector<Real> fixed_occupancy_;
    std::vector<Real> occupancy_;
    std::vector<Real> excess_density_;

    void add_rectangle(std::vector<Real>& map, const Node& node) const;
    void append_rectangle_changes(const Node& node, Real x, Real y, Real sign,
                                  std::vector<std::pair<int, Real>>& changes) const;
};

}  // namespace ea
