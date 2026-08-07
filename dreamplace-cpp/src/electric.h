#pragma once

#include "types.h"

#include <cstdint>
#include <vector>

namespace dpcpp {

struct Filler {
    Real x = 0.0;
    Real y = 0.0;
    Real width = 0.0;
    Real height = 0.0;
};

struct DensityResult {
    Real energy = 0.0;
    Real overflow = 0.0;
    Real max_density = 0.0;
};

std::vector<Filler> initialize_fillers(const Database& db, Real target_density,
                                       std::uint64_t seed);

class ElectricDensity {
public:
    ElectricDensity(const Database& db, int bins_x, int bins_y,
                    Real target_density);

    DensityResult compute(const Database& db, const std::vector<Filler>& fillers,
                          std::vector<Real>* node_grad_x,
                          std::vector<Real>* node_grad_y,
                          std::vector<Real>* filler_grad_x,
                          std::vector<Real>* filler_grad_y);

    int bins_x() const { return nx_; }
    int bins_y() const { return ny_; }
    Real bin_width() const { return bin_w_; }
    Real bin_height() const { return bin_h_; }
    Real target_density() const { return target_density_; }
    const std::vector<Real>& density_map() const { return density_; }

private:
    void deposit_rectangle(Real cx, Real cy, Real width, Real height,
                           Real ratio, std::vector<Real>& map) const;
    void gather_gradient(Real cx, Real cy, Real width, Real height, Real ratio,
                         const std::vector<Real>& grad_phi_x,
                         const std::vector<Real>& grad_phi_y,
                         Real& gx, Real& gy) const;

    int nx_ = 0;
    int ny_ = 0;
    Real xl_ = 0.0;
    Real yl_ = 0.0;
    Real xh_ = 0.0;
    Real yh_ = 0.0;
    Real bin_w_ = 0.0;
    Real bin_h_ = 0.0;
    Real target_density_ = 1.0;
    Real movable_area_ = 0.0;
    std::vector<Real> fixed_density_;
    std::vector<Real> density_;
};

void clamp_to_region(Database& db, std::vector<Filler>& fillers);

}  // namespace dpcpp

