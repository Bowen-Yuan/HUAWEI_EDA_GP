#include "../src/types.h"

#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

void density_init(BinGrid&, float_t, float_t, float_t, float_t,
                  const std::vector<Cell>&, float_t);
float_t density_compute_and_gradient(const BinGrid&, const std::vector<Cell>&,
                                     std::vector<float_t>&, std::vector<float_t>&,
                                     std::vector<float_t>&);

int main() {
    std::vector<Cell> cells = {
        {0, 7.5f, 8.0f, 1.0f, 1.0f, false, 1.0f},
        {1, 8.5f, 8.0f, 1.0f, 1.0f, false, 1.0f},
    };
    BinGrid grid{};
    density_init(grid, 0.0f, 0.0f, 16.0f, 16.0f, cells, 2.0f / 256.0f);

    std::vector<float_t> gx(cells.size(), 0.0f);
    std::vector<float_t> gy(cells.size(), 0.0f);
    std::vector<float_t> density;
    const float_t energy = density_compute_and_gradient(grid, cells, gx, gy, density);
    const float_t deposited_area = std::accumulate(density.begin(), density.end(), 0.0f);

    if (!std::isfinite(energy) || energy <= 0.0f) {
        std::fprintf(stderr, "invalid electric energy: %g\n", energy);
        return 1;
    }
    if (std::fabs(deposited_area - 2.0f) > 1e-4f) {
        std::fprintf(stderr, "charge is not conserved: %g\n", deposited_area);
        return 2;
    }
    if (!(gx[0] > 0.0f && gx[1] < 0.0f)) {
        std::fprintf(stderr, "crowded cells are not pushed apart: gx=(%g,%g)\n",
                     gx[0], gx[1]);
        return 3;
    }
    if (std::fabs(gx[0]) < 1e-6f || std::fabs(gx[1]) < 1e-6f) {
        std::fprintf(stderr, "electric gradient unexpectedly vanished\n");
        return 4;
    }

    std::printf("electric density test passed: energy=%g area=%g gx=(%g,%g)\n",
                energy, deposited_area, gx[0], gx[1]);
    return 0;
}
