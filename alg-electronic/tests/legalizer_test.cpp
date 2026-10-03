#include "../src/legalizer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

namespace {

bool overlaps(double a0, double a1, double b0, double b1) {
    return a0 < b1 - 1.0e-5 && b0 < a1 - 1.0e-5;
}

void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "Legalizer test failed: %s\n", message);
        std::exit(1);
    }
}

}  // namespace

int main() {
    const char* scl_file = "build/legalizer_test.scl";
    {
        std::ofstream out(scl_file);
        out << "UCLA scl 1.0\nNumRows : 2\n";
        for (int y : {0, 10}) {
            out << "CoreRow Horizontal\n"
                << "  Coordinate : " << y << "\n"
                << "  Height : 10\n"
                << "  Sitewidth : 1\n"
                << "  Sitespacing : 1\n"
                << "  Siteorient : N\n"
                << "  Sitesymmetry : Y\n"
                << "  SubrowOrigin : 0 NumSites : 30\n"
                << "End\n";
        }
    }

    std::vector<Cell> cells;
    for (int i = 0; i < 8; ++i) {
        cells.push_back(Cell{i, 15.0f, 10.0f, 5.0f, 10.0f, false, 50.0f});
    }
    // Fixed macro cuts the interval [12,18) out of both rows.
    cells.push_back(Cell{8, 15.0f, 10.0f, 6.0f, 20.0f, true, 120.0f});

    const std::vector<Net> nets;
    LegalizerConfig config;
    config.mode = LegalizerMode::DreamplaceGreedyAbacus;
    config.fallback_to_legacy = false;
    const auto result = legalize_placement(cells, nets, 0.0f, 0.0f, 30.0f,
                                           20.0f, scl_file, config);
    require(result.success, result.message.c_str());
    require(result.movable_cells == 8, "movable-cell count mismatch");
    require(cells[8].x == 15.0f && cells[8].y == 10.0f,
            "fixed macro moved");

    for (int i = 0; i < 8; ++i) {
        const double xl = cells[i].x - 0.5 * cells[i].width;
        const double xh = cells[i].x + 0.5 * cells[i].width;
        const double yl = cells[i].y - 0.5 * cells[i].height;
        require(std::abs(xl - std::round(xl)) < 1.0e-5,
                "cell is not site-aligned");
        require(std::abs(yl) < 1.0e-5 || std::abs(yl - 10.0) < 1.0e-5,
                "cell is not row-aligned");
        require(xl >= -1.0e-5 && xh <= 30.0 + 1.0e-5,
                "cell is outside the row boundary");
        require(!overlaps(xl, xh, 12.0, 18.0),
                "cell overlaps the fixed macro");
        for (int j = 0; j < i; ++j) {
            const double jxl = cells[j].x - 0.5 * cells[j].width;
            const double jxh = cells[j].x + 0.5 * cells[j].width;
            const double jyl = cells[j].y - 0.5 * cells[j].height;
            require(std::abs(yl - jyl) > 1.0e-5 ||
                    !overlaps(xl, xh, jxl, jxh),
                    "movable cells overlap");
        }
    }

    std::printf("Legalizer test passed: %zu cells, %zu free segments\n",
                result.movable_cells, result.free_segments);

    // The original public call remains available and keeps its legacy behavior.
    for (int i = 0; i < 8; ++i) {
        cells[i].x = 15.0f;
        cells[i].y = 10.0f;
    }
    const auto legacy = legalize_placement(cells, nets, 0.0f, 0.0f, 30.0f,
                                            20.0f, scl_file);
    require(legacy.success, "legacy compatibility path failed");
    return 0;
}
