#include "types.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include <omp.h>

namespace {

constexpr float_t kSqrt2 = 1.4142135623730951f;
constexpr double kPi = 3.1415926535897932384626433832795;

// Kept for the legacy spread-only solver. The regular placement solvers use
// an unscaled electric gradient and apply lambda explicitly at the call site.
float_t g_density_grad_scale = 1.0f;
bool g_short_edge_only = false;

struct AxisWeights {
    int begin = 0;
    std::vector<float_t> values;
};

struct SpectralCache {
    int nx = 0;
    int ny = 0;
    float_t width = 0.0f;
    float_t height = 0.0f;
    std::vector<double> cos_x;
    std::vector<double> sin_x;
    std::vector<double> cos_y;
    std::vector<double> sin_y;
    std::vector<double> kx;
    std::vector<double> ky;
};

std::vector<SpectralCache> g_spectral_caches;

AxisWeights make_axis_weights(float_t center, float_t node_size,
                              float_t lower, float_t bin_size, int num_bins) {
    AxisWeights result;
    if (num_bins <= 0 || bin_size <= 0.0f) return result;

    // DREAMPlace stretches small cells so that the density function remains
    // smooth on the bin grid. A one-bin linear shoulder is then used as a
    // compact triangular smoothing kernel around the stretched rectangle.
    const float_t clamped_size = std::max(node_size, kSqrt2 * bin_size);
    const float_t half_size = 0.5f * clamped_size;
    const float_t shoulder = kSqrt2 * bin_size;
    const float_t support = half_size + shoulder;

    int first = (int)std::floor((center - support - lower) / bin_size);
    int last = (int)std::floor((center + support - lower) / bin_size);
    first = std::max(0, first);
    last = std::min(num_bins - 1, last);

    if (first > last) {
        const int nearest = std::max(0, std::min(num_bins - 1,
            (int)std::floor((center - lower) / bin_size)));
        result.begin = nearest;
        result.values.assign(1, 1.0f);
        return result;
    }

    result.begin = first;
    result.values.resize(last - first + 1, 0.0f);
    double sum = 0.0;
    for (int b = first; b <= last; ++b) {
        const float_t bin_center = lower + (b + 0.5f) * bin_size;
        const float_t outside = std::max(0.0f, std::fabs(bin_center - center) - half_size);
        const float_t weight = std::max(0.0f, 1.0f - outside / shoulder);
        result.values[b - first] = weight;
        sum += weight;
    }

    if (sum <= 0.0) {
        const int nearest = std::max(first, std::min(last,
            (int)std::floor((center - lower) / bin_size)));
        std::fill(result.values.begin(), result.values.end(), 0.0f);
        result.values[nearest - first] = 1.0f;
    } else {
        const float_t inv_sum = (float_t)(1.0 / sum);
        for (float_t& weight : result.values) weight *= inv_sum;
    }
    return result;
}

void deposit_movable_cell(const BinGrid& grid, const Cell& cell,
                          std::vector<float_t>& density) {
    const AxisWeights wx = make_axis_weights(cell.x, cell.width,
        grid.chip_xl, grid.bin_w, grid.nx);
    const AxisWeights wy = make_axis_weights(cell.y, cell.height,
        grid.chip_yl, grid.bin_h, grid.ny);

    for (int iy = 0; iy < (int)wy.values.size(); ++iy) {
        const int by = wy.begin + iy;
        for (int ix = 0; ix < (int)wx.values.size(); ++ix) {
            const int bx = wx.begin + ix;
            density[by * grid.nx + bx] +=
                cell.area * wx.values[ix] * wy.values[iy];
        }
    }
}

void deposit_fixed_cell(const BinGrid& grid, const Cell& cell,
                        std::vector<float_t>& density) {
    const float_t cl = cell.x - 0.5f * cell.width;
    const float_t cr = cell.x + 0.5f * cell.width;
    const float_t cb = cell.y - 0.5f * cell.height;
    const float_t ct = cell.y + 0.5f * cell.height;

    int bx0 = std::max(0, (int)std::floor((cl - grid.chip_xl) / grid.bin_w));
    int bx1 = std::min(grid.nx - 1, (int)std::floor((cr - grid.chip_xl) / grid.bin_w));
    int by0 = std::max(0, (int)std::floor((cb - grid.chip_yl) / grid.bin_h));
    int by1 = std::min(grid.ny - 1, (int)std::floor((ct - grid.chip_yl) / grid.bin_h));
    if (bx0 > bx1 || by0 > by1) return;

    for (int by = by0; by <= by1; ++by) {
        const float_t byl = grid.chip_yl + by * grid.bin_h;
        const float_t byh = byl + grid.bin_h;
        const float_t oy = std::max(0.0f, std::min(ct, byh) - std::max(cb, byl));
        for (int bx = bx0; bx <= bx1; ++bx) {
            const float_t bxl = grid.chip_xl + bx * grid.bin_w;
            const float_t bxh = bxl + grid.bin_w;
            const float_t ox = std::max(0.0f, std::min(cr, bxh) - std::max(cl, bxl));
            // Match DREAMPlace: fixed blockage density is scaled by target density.
            density[by * grid.nx + bx] += grid.target_density * ox * oy;
        }
    }
}

void build_density_map(const BinGrid& grid, const std::vector<Cell>& cells,
                       std::vector<float_t>& density) {
    const int nbins = grid.nx * grid.ny;
    density.assign(nbins, 0.0f);

    #pragma omp parallel
    {
        std::vector<float_t> local(nbins, 0.0f);
        #pragma omp for schedule(dynamic, 64)
        for (int i = 0; i < (int)cells.size(); ++i) {
            if (cells[i].is_terminal)
                deposit_fixed_cell(grid, cells[i], local);
            else
                deposit_movable_cell(grid, cells[i], local);
        }
        #pragma omp critical
        {
            for (int b = 0; b < nbins; ++b) density[b] += local[b];
        }
    }
}

SpectralCache& spectral_cache(const BinGrid& grid) {
    const float_t width = grid.chip_xh - grid.chip_xl;
    const float_t height = grid.chip_yh - grid.chip_yl;
    for (SpectralCache& cache : g_spectral_caches) {
        if (cache.nx == grid.nx && cache.ny == grid.ny &&
            cache.width == width && cache.height == height) return cache;
    }

    SpectralCache cache;
    cache.nx = grid.nx;
    cache.ny = grid.ny;
    cache.width = width;
    cache.height = height;
    cache.cos_x.resize(grid.nx * grid.nx);
    cache.sin_x.resize(grid.nx * grid.nx);
    cache.cos_y.resize(grid.ny * grid.ny);
    cache.sin_y.resize(grid.ny * grid.ny);
    cache.kx.resize(grid.nx);
    cache.ky.resize(grid.ny);

    for (int u = 0; u < grid.nx; ++u) {
        const double alpha = (u == 0) ? std::sqrt(1.0 / grid.nx)
                                      : std::sqrt(2.0 / grid.nx);
        cache.kx[u] = kPi * u / width;
        for (int x = 0; x < grid.nx; ++x) {
            const double angle = kPi * u * (x + 0.5) / grid.nx;
            cache.cos_x[u * grid.nx + x] = alpha * std::cos(angle);
            cache.sin_x[u * grid.nx + x] = alpha * std::sin(angle);
        }
    }
    for (int v = 0; v < grid.ny; ++v) {
        const double alpha = (v == 0) ? std::sqrt(1.0 / grid.ny)
                                      : std::sqrt(2.0 / grid.ny);
        cache.ky[v] = kPi * v / height;
        for (int y = 0; y < grid.ny; ++y) {
            const double angle = kPi * v * (y + 0.5) / grid.ny;
            cache.cos_y[v * grid.ny + y] = alpha * std::cos(angle);
            cache.sin_y[v * grid.ny + y] = alpha * std::sin(angle);
        }
    }

    g_spectral_caches.push_back(std::move(cache));
    return g_spectral_caches.back();
}

double solve_electric_field(const BinGrid& grid,
                            const std::vector<float_t>& bin_density,
                            std::vector<float_t>& field_x,
                            std::vector<float_t>& field_y) {
    const int nx = grid.nx;
    const int ny = grid.ny;
    const int nbins = nx * ny;
    const double bin_area = (double)grid.bin_w * grid.bin_h;
    SpectralCache& cache = spectral_cache(grid);

    std::vector<double> rho(nbins, 0.0);
    for (int b = 0; b < nbins; ++b) rho[b] = bin_density[b] / bin_area;

    // Orthonormal 2D DCT-II: rho(x,y) -> rho_hat(u,v).
    std::vector<double> temp(nx * ny, 0.0);
    std::vector<double> rho_hat(nbins, 0.0);
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; ++y) {
        for (int u = 0; u < nx; ++u) {
            double sum = 0.0;
            const double* basis = &cache.cos_x[u * nx];
            for (int x = 0; x < nx; ++x) sum += basis[x] * rho[y * nx + x];
            temp[y * nx + u] = sum;
        }
    }
    #pragma omp parallel for schedule(static)
    for (int v = 0; v < ny; ++v) {
        const double* basis = &cache.cos_y[v * ny];
        for (int u = 0; u < nx; ++u) {
            double sum = 0.0;
            for (int y = 0; y < ny; ++y) sum += basis[y] * temp[y * nx + u];
            rho_hat[v * nx + u] = sum;
        }
    }

    // Spectral Poisson solve. The DC component is removed, which is
    // equivalent to adding the uniform neutralizing background charge.
    std::vector<double> phi_hat(nbins, 0.0);
    double energy = 0.0;
    for (int v = 0; v < ny; ++v) {
        for (int u = 0; u < nx; ++u) {
            if (u == 0 && v == 0) continue;
            const double k2 = cache.kx[u] * cache.kx[u] +
                              cache.ky[v] * cache.ky[v];
            const int index = v * nx + u;
            phi_hat[index] = rho_hat[index] / k2;
            energy += 0.5 * bin_area * rho_hat[index] * phi_hat[index];
        }
    }

    // E_x = -d(phi)/dx and E_y = -d(phi)/dy. Differentiating the cosine
    // basis produces the mixed sine/cosine inverse transforms below.
    std::vector<double> mixed(nx * ny, 0.0);
    field_x.assign(nbins, 0.0f);
    field_y.assign(nbins, 0.0f);

    #pragma omp parallel for schedule(static)
    for (int v = 0; v < ny; ++v) {
        for (int x = 0; x < nx; ++x) {
            double sum = 0.0;
            for (int u = 1; u < nx; ++u) {
                sum += cache.sin_x[u * nx + x] * cache.kx[u] *
                       phi_hat[v * nx + u];
            }
            mixed[v * nx + x] = sum;
        }
    }
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            double sum = 0.0;
            for (int v = 0; v < ny; ++v)
                sum += cache.cos_y[v * ny + y] * mixed[v * nx + x];
            field_x[y * nx + x] = (float_t)sum;
        }
    }

    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; ++y) {
        for (int u = 0; u < nx; ++u) {
            double sum = 0.0;
            for (int v = 1; v < ny; ++v) {
                sum += cache.sin_y[v * ny + y] * cache.ky[v] *
                       phi_hat[v * nx + u];
            }
            temp[y * nx + u] = sum;
        }
    }
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            double sum = 0.0;
            for (int u = 0; u < nx; ++u)
                sum += cache.cos_x[u * nx + x] * temp[y * nx + u];
            field_y[y * nx + x] = (float_t)sum;
        }
    }
    return energy;
}

void accumulate_cell_gradient(const BinGrid& grid, const Cell& cell,
                              const std::vector<float_t>& field_x,
                              const std::vector<float_t>& field_y,
                              float_t& grad_x, float_t& grad_y) {
    const AxisWeights wx = make_axis_weights(cell.x, cell.width,
        grid.chip_xl, grid.bin_w, grid.nx);
    const AxisWeights wy = make_axis_weights(cell.y, cell.height,
        grid.chip_yl, grid.bin_h, grid.ny);
    double force_x = 0.0;
    double force_y = 0.0;
    for (int iy = 0; iy < (int)wy.values.size(); ++iy) {
        const int by = wy.begin + iy;
        for (int ix = 0; ix < (int)wx.values.size(); ++ix) {
            const int bx = wx.begin + ix;
            const double area_weight = (double)cell.area * wx.values[ix] * wy.values[iy];
            const int index = by * grid.nx + bx;
            force_x += area_weight * field_x[index];
            force_y += area_weight * field_y[index];
        }
    }
    // The density-energy gradient is negative electric force. Gradient descent
    // therefore moves positive charges along the field, out of crowded areas.
    float_t density_grad_x = -g_density_grad_scale * (float_t)force_x;
    float_t density_grad_y = -g_density_grad_scale * (float_t)force_y;
    if (g_short_edge_only) {
        if (std::fabs(density_grad_x) > std::fabs(density_grad_y))
            density_grad_y = 0.0f;
        else
            density_grad_x = 0.0f;
    }
    grad_x += density_grad_x;
    grad_y += density_grad_y;
}

}  // namespace

void set_density_grad_scale(float_t scale) { g_density_grad_scale = scale; }
float_t get_density_grad_scale() { return g_density_grad_scale; }
void set_short_edge_only(bool enabled) { g_short_edge_only = enabled; }

void density_init(BinGrid& grid,
                  float_t chip_xl, float_t chip_yl,
                  float_t chip_xh, float_t chip_yh,
                  const std::vector<Cell>& cells,
                  float_t target_density) {
    grid.chip_xl = chip_xl;
    grid.chip_yl = chip_yl;
    grid.chip_xh = chip_xh;
    grid.chip_yh = chip_yh;
    grid.target_density = target_density;

    const float_t chip_w = chip_xh - chip_xl;
    const float_t chip_h = chip_yh - chip_yl;
    double fixed_blockage_area = 0.0;
    for (const Cell& cell : cells) {
        if (!cell.is_terminal) continue;
        const double left = std::max((double)chip_xl,
                                     (double)cell.x - 0.5 * cell.width);
        const double right = std::min((double)chip_xh,
                                      (double)cell.x + 0.5 * cell.width);
        const double bottom = std::max((double)chip_yl,
                                       (double)cell.y - 0.5 * cell.height);
        const double top = std::min((double)chip_yh,
                                    (double)cell.y + 0.5 * cell.height);
        if (right > left && top > bottom)
            fixed_blockage_area += (right - left) * (top - bottom);
    }
    const double chip_area = (double)chip_w * chip_h;
    grid.available_area = (float_t)std::max(chip_area * 0.01,
                                             chip_area - fixed_blockage_area);
    float_t avg_w = 0.0f;
    float_t avg_h = 0.0f;
    int movable_count = 0;
    for (const Cell& cell : cells) {
        if (!cell.is_terminal) {
            avg_w += cell.width;
            avg_h += cell.height;
            ++movable_count;
        }
    }
    if (movable_count == 0) {
        avg_w = chip_w / 100.0f;
        avg_h = chip_h / 100.0f;
    } else {
        avg_w /= movable_count;
        avg_h /= movable_count;
    }

    const float_t target_bin_w = std::max(avg_w * 15.0f, chip_w / 200.0f);
    const float_t target_bin_h = std::max(avg_h * 15.0f, chip_h / 200.0f);
    grid.nx = std::max(8, (int)(chip_w / target_bin_w));
    grid.ny = std::max(8, (int)(chip_h / target_bin_h));
    grid.bin_w = chip_w / grid.nx;
    grid.bin_h = chip_h / grid.ny;

    const int nbins = grid.nx * grid.ny;
    grid.density.assign(nbins, 0.0f);
    grid.overflow.assign(nbins, 0.0f);
    grid.total_overflow = 0.0f;
    std::printf("[ElectricDensity] Grid: %d x %d = %d bins, bin size %.1f x %.1f, target=%.4f\n",
                grid.nx, grid.ny, nbins, grid.bin_w, grid.bin_h, target_density);
}

void density_build_map(const BinGrid& grid, const std::vector<Cell>& cells,
                       std::vector<float_t>& density) {
    build_density_map(grid, cells, density);
}

float_t density_compute_and_gradient(
    const BinGrid& grid,
    const std::vector<Cell>& cells,
    std::vector<float_t>& grad_x,
    std::vector<float_t>& grad_y,
    std::vector<float_t>& bin_density) {
    build_density_map(grid, cells, bin_density);

    std::vector<float_t> field_x;
    std::vector<float_t> field_y;
    const double energy = solve_electric_field(grid, bin_density, field_x, field_y);

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < (int)cells.size(); ++i) {
        if (cells[i].is_terminal) continue;
        accumulate_cell_gradient(grid, cells[i], field_x, field_y,
                                 grad_x[i], grad_y[i]);
    }
    return (float_t)energy;
}

float_t density_evaluate(const BinGrid& grid, const std::vector<Cell>& cells) {
    std::vector<float_t> bins;
    build_density_map(grid, cells, bins);
    const float_t inv_bin_area = 1.0f / (grid.bin_w * grid.bin_h);
    double excess_area = 0.0;
    float_t max_overflow = 0.0f;
    for (float_t area : bins) {
        const float_t overflow = area * inv_bin_area - grid.target_density;
        if (overflow > 0.0f) {
            excess_area += (double)overflow * grid.bin_w * grid.bin_h;
            max_overflow = std::max(max_overflow, overflow);
        }
    }
    const float_t overflow_ratio = grid.available_area > 0.0f
        ? (float_t)(excess_area / grid.available_area) : 0.0f;
    std::printf("[ElectricDensity] overflow_ratio=%.6f, max_overflow=%.6f\n",
                overflow_ratio, max_overflow);
    return overflow_ratio;
}
