#include "electric.h"

#include "spectral.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

namespace dpcpp {
namespace {
constexpr Real kPi = 3.141592653589793238462643383279502884;

Real free_row_area(const Database& db) {
    Real total = 0.0;
    for (const Row& row : db.rows) {
        std::vector<std::pair<Real, Real>> blocked;
        for (int id : db.fixed_ids) {
            const Node& node = db.nodes[id];
            if (node.terminal_ni) continue;
            const Real bottom = node.y - 0.5 * node.height;
            const Real top = node.y + 0.5 * node.height;
            if (top <= row.y || bottom >= row.y + row.height) continue;
            const Real left = std::max(row.origin, node.x - 0.5 * node.width);
            const Real right = std::min(row.xh(), node.x + 0.5 * node.width);
            if (right > left) blocked.emplace_back(left, right);
        }
        std::sort(blocked.begin(), blocked.end());
        Real blocked_width = 0.0;
        Real begin = 0.0, end = 0.0;
        bool active = false;
        for (const auto& interval : blocked) {
            if (!active || interval.first > end) {
                if (active) blocked_width += end - begin;
                begin = interval.first;
                end = interval.second;
                active = true;
            } else {
                end = std::max(end, interval.second);
            }
        }
        if (active) blocked_width += end - begin;
        total += std::max<Real>(0.0, row.xh() - row.origin - blocked_width) * row.height;
    }
    return total;
}

Real standard_cell_width(const Database& db) {
    const Real row_height = db.rows.front().height;
    std::vector<Real> widths;
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        if (node.height <= row_height * 1.5) widths.push_back(node.width);
    }
    if (widths.empty()) return db.rows.front().site_spacing;
    std::sort(widths.begin(), widths.end());
    const std::size_t lo = widths.size() / 20;
    const std::size_t hi = widths.size() - widths.size() / 20;
    Real sum = 0.0;
    for (std::size_t i = lo; i < hi; ++i) sum += widths[i];
    return sum / std::max<std::size_t>(1, hi - lo);
}

}  // namespace

std::vector<Filler> initialize_fillers(const Database& db, Real target_density,
                                       std::uint64_t seed) {
    Real movable_area = 0.0;
    for (int id : db.movable_ids) movable_area += db.nodes[id].area();
    const Real placeable = free_row_area(db);
    const Real filler_area_total = std::max<Real>(0.0, placeable * target_density - movable_area);
    const Real filler_w = standard_cell_width(db);
    const Real filler_h = db.rows.front().height;
    const Real one_area = filler_w * filler_h;
    const std::size_t count = one_area > 0.0
        ? static_cast<std::size_t>(std::llround(filler_area_total / one_area)) : 0;
    std::mt19937_64 generator(seed ^ 0x9e3779b97f4a7c15ULL);
    std::uniform_real_distribution<Real> ux(db.xl + 0.5 * filler_w,
                                            db.xh - 0.5 * filler_w);
    std::uniform_real_distribution<Real> uy(db.yl + 0.5 * filler_h,
                                            db.yh - 0.5 * filler_h);
    std::vector<Filler> fillers(count);
    for (Filler& filler : fillers) {
        filler.x = ux(generator);
        filler.y = uy(generator);
        filler.width = filler_w;
        filler.height = filler_h;
    }
    std::cout << "[Fillers] free_row_area=" << placeable
              << " movable_area=" << movable_area
              << " filler_area=" << filler_area_total
              << " count=" << count
              << " size=" << filler_w << 'x' << filler_h << '\n';
    return fillers;
}

ElectricDensity::ElectricDensity(const Database& db, int bins_x, int bins_y,
                                  Real target_density, ElectricFieldModel field_model)
    : nx_(bins_x), ny_(bins_y), xl_(db.xl), yl_(db.yl), xh_(db.xh), yh_(db.yh),
      target_density_(target_density), field_model_(field_model) {
    if (!is_power_of_two(nx_) || !is_power_of_two(ny_)) {
        throw std::runtime_error("electric grid dimensions must be powers of two");
    }
    bin_w_ = (xh_ - xl_) / nx_;
    bin_h_ = (yh_ - yl_) / ny_;
    fixed_density_.assign(nx_ * ny_, 0.0);
    for (int id : db.fixed_ids) {
        const Node& node = db.nodes[id];
        if (node.terminal_ni) continue;
        deposit_rectangle(node.x, node.y, node.width, node.height,
                          target_density_, fixed_density_);
    }
    for (int id : db.movable_ids) movable_area_ += db.nodes[id].area();
    std::cout << "[Electric] grid=" << nx_ << 'x' << ny_
              << " bin=" << bin_w_ << 'x' << bin_h_
              << " target=" << target_density_
              << " field=" << (field_model_ == ElectricFieldModel::MixedSpectral
                                    ? "mixed-spectral" : "finite-difference") << '\n';
}

void ElectricDensity::deposit_rectangle(Real cx, Real cy, Real width, Real height,
                                        Real ratio, std::vector<Real>& map) const {
    const Real left = cx - 0.5 * width;
    const Real right = cx + 0.5 * width;
    const Real bottom = cy - 0.5 * height;
    const Real top = cy + 0.5 * height;
    const int bx0 = std::max(0, static_cast<int>(std::floor((left - xl_) / bin_w_)));
    const int bx1 = std::min(nx_ - 1, static_cast<int>(std::floor((right - xl_) / bin_w_)));
    const int by0 = std::max(0, static_cast<int>(std::floor((bottom - yl_) / bin_h_)));
    const int by1 = std::min(ny_ - 1, static_cast<int>(std::floor((top - yl_) / bin_h_)));
    for (int by = by0; by <= by1; ++by) {
        const Real bin_bottom = yl_ + by * bin_h_;
        const Real overlap_y = std::max<Real>(0.0,
            std::min(top, bin_bottom + bin_h_) - std::max(bottom, bin_bottom));
        for (int bx = bx0; bx <= bx1; ++bx) {
            const Real bin_left = xl_ + bx * bin_w_;
            const Real overlap_x = std::max<Real>(0.0,
                std::min(right, bin_left + bin_w_) - std::max(left, bin_left));
            map[by * nx_ + bx] += ratio * overlap_x * overlap_y;
        }
    }
}

void ElectricDensity::gather_gradient(Real cx, Real cy, Real width, Real height,
                                      Real ratio,
                                      const std::vector<Real>& grad_phi_x,
                                      const std::vector<Real>& grad_phi_y,
                                      Real& gx, Real& gy) const {
    const Real left = cx - 0.5 * width;
    const Real right = cx + 0.5 * width;
    const Real bottom = cy - 0.5 * height;
    const Real top = cy + 0.5 * height;
    const int bx0 = std::max(0, static_cast<int>(std::floor((left - xl_) / bin_w_)));
    const int bx1 = std::min(nx_ - 1, static_cast<int>(std::floor((right - xl_) / bin_w_)));
    const int by0 = std::max(0, static_cast<int>(std::floor((bottom - yl_) / bin_h_)));
    const int by1 = std::min(ny_ - 1, static_cast<int>(std::floor((top - yl_) / bin_h_)));
    gx = 0.0;
    gy = 0.0;
    for (int by = by0; by <= by1; ++by) {
        const Real bin_bottom = yl_ + by * bin_h_;
        const Real overlap_y = std::max<Real>(0.0,
            std::min(top, bin_bottom + bin_h_) - std::max(bottom, bin_bottom));
        for (int bx = bx0; bx <= bx1; ++bx) {
            const Real bin_left = xl_ + bx * bin_w_;
            const Real overlap_x = std::max<Real>(0.0,
                std::min(right, bin_left + bin_w_) - std::max(left, bin_left));
            const Real charge = ratio * overlap_x * overlap_y;
            gx += charge * grad_phi_x[by * nx_ + bx];
            gy += charge * grad_phi_y[by * nx_ + bx];
        }
    }
}

DensityResult ElectricDensity::compute(
    const Database& db, const std::vector<Filler>& fillers,
    std::vector<Real>* node_grad_x, std::vector<Real>* node_grad_y,
    std::vector<Real>* filler_grad_x, std::vector<Real>* filler_grad_y) {
    density_ = fixed_density_;
    const Real sqrt2 = std::sqrt(2.0);
    for (int id : db.movable_ids) {
        const Node& node = db.nodes[id];
        const Real width = std::max(node.width, sqrt2 * bin_w_);
        const Real height = std::max(node.height, sqrt2 * bin_h_);
        const Real ratio = node.area() / (width * height);
        deposit_rectangle(node.x, node.y, width, height, ratio, density_);
    }
    const Real bin_area = bin_w_ * bin_h_;
    DensityResult result;
    // DREAMPlace evaluates overflow with num_filler_nodes=0.  Fillers belong
    // to the electrostatic objective, but including them in the stopping
    // metric biases overflow upward and corrupts gamma/checkpoint selection.
    for (int i = 0; i < nx_ * ny_; ++i) {
        result.max_density = std::max(result.max_density, density_[i] / bin_area);
        result.overflow += std::max<Real>(
            0.0, density_[i] - target_density_ * bin_area);
    }
    result.overflow /= std::max<Real>(1.0, movable_area_);

    for (const Filler& filler : fillers) {
        const Real width = std::max(filler.width, sqrt2 * bin_w_);
        const Real height = std::max(filler.height, sqrt2 * bin_h_);
        const Real ratio = filler.width * filler.height / (width * height);
        deposit_rectangle(filler.x, filler.y, width, height, ratio, density_);
    }

    std::vector<Real> rho(nx_ * ny_);
    for (int i = 0; i < nx_ * ny_; ++i) {
        const Real relative = density_[i] / bin_area;
        rho[i] = relative;
    }

    std::vector<Real> coefficients = rho;
    dct2_orthonormal(coefficients, nx_, ny_);
    for (int y = 0; y < ny_; ++y) {
        const Real ky = kPi * y / (yh_ - yl_);
        for (int x = 0; x < nx_; ++x) {
            const Real kx = kPi * x / (xh_ - xl_);
            const int index = y * nx_ + x;
            const Real k2 = kx * kx + ky * ky;
            coefficients[index] = k2 > 0.0 ? coefficients[index] / k2 : 0.0;
        }
    }
    const std::vector<Real> potential_coefficients = coefficients;
    std::vector<Real> potential = potential_coefficients;
    idct2_orthonormal(potential, nx_, ny_);
    for (int i = 0; i < nx_ * ny_; ++i) result.energy += 0.5 * rho[i] * potential[i] * bin_area;

    if (!node_grad_x && !filler_grad_x) return result;
    std::vector<Real> grad_phi_x(nx_ * ny_, 0.0);
    std::vector<Real> grad_phi_y(nx_ * ny_, 0.0);
    if (field_model_ == ElectricFieldModel::MixedSpectral) {
        grad_phi_x = potential_coefficients;
        grad_phi_y = potential_coefficients;
        #pragma omp parallel for schedule(static)
        for (int y = 0; y < ny_; ++y) {
            const Real ky = kPi * y / (yh_ - yl_);
            for (int x = 0; x < nx_; ++x) {
                const Real kx = kPi * x / (xh_ - xl_);
                const int index = y * nx_ + x;
                grad_phi_x[index] *= -kx;
                grad_phi_y[index] *= -ky;
            }
        }
        inverse_mixed_sine_cosine2(grad_phi_x, nx_, ny_, true);
        inverse_mixed_sine_cosine2(grad_phi_y, nx_, ny_, false);
    } else {
        #pragma omp parallel for schedule(static)
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const int xm = std::max(0, x - 1);
                const int xp = std::min(nx_ - 1, x + 1);
                const int ym = std::max(0, y - 1);
                const int yp = std::min(ny_ - 1, y + 1);
                grad_phi_x[y * nx_ + x] =
                    (potential[y * nx_ + xp] - potential[y * nx_ + xm]) /
                    ((xp - xm) * bin_w_);
                grad_phi_y[y * nx_ + x] =
                    (potential[yp * nx_ + x] - potential[ym * nx_ + x]) /
                    ((yp - ym) * bin_h_);
            }
        }
    }
    if (node_grad_x) node_grad_x->assign(db.nodes.size(), 0.0);
    if (node_grad_y) node_grad_y->assign(db.nodes.size(), 0.0);
    #pragma omp parallel for schedule(dynamic, 256)
    for (int m = 0; m < static_cast<int>(db.movable_ids.size()); ++m) {
        const int id = db.movable_ids[m];
        const Node& node = db.nodes[id];
        const Real width = std::max(node.width, sqrt2 * bin_w_);
        const Real height = std::max(node.height, sqrt2 * bin_h_);
        const Real ratio = node.area() / (width * height);
        Real gx = 0.0, gy = 0.0;
        gather_gradient(node.x, node.y, width, height, ratio,
                        grad_phi_x, grad_phi_y, gx, gy);
        if (node_grad_x) (*node_grad_x)[id] = gx;
        if (node_grad_y) (*node_grad_y)[id] = gy;
    }
    if (filler_grad_x) filler_grad_x->assign(fillers.size(), 0.0);
    if (filler_grad_y) filler_grad_y->assign(fillers.size(), 0.0);
    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < static_cast<int>(fillers.size()); ++i) {
        const Filler& filler = fillers[i];
        const Real width = std::max(filler.width, sqrt2 * bin_w_);
        const Real height = std::max(filler.height, sqrt2 * bin_h_);
        const Real ratio = filler.width * filler.height / (width * height);
        Real gx = 0.0, gy = 0.0;
        gather_gradient(filler.x, filler.y, width, height, ratio,
                        grad_phi_x, grad_phi_y, gx, gy);
        if (filler_grad_x) (*filler_grad_x)[i] = gx;
        if (filler_grad_y) (*filler_grad_y)[i] = gy;
    }
    return result;
}

void clamp_to_region(Database& db, std::vector<Filler>& fillers) {
    for (int id : db.movable_ids) {
        Node& node = db.nodes[id];
        node.x = std::clamp(node.x, db.xl + 0.5 * node.width,
                            db.xh - 0.5 * node.width);
        node.y = std::clamp(node.y, db.yl + 0.5 * node.height,
                            db.yh - 0.5 * node.height);
    }
    for (Filler& filler : fillers) {
        filler.x = std::clamp(filler.x, db.xl + 0.5 * filler.width,
                              db.xh - 0.5 * filler.width);
        filler.y = std::clamp(filler.y, db.yl + 0.5 * filler.height,
                              db.yh - 0.5 * filler.height);
    }
}

}  // namespace dpcpp
