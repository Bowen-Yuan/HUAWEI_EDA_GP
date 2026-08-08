#include "bookshelf.h"
#include "legalizer.h"
#include "spectral.h"
#include "wirelength.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace dpcpp;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    {
        Database db;
        db.xl = 0.0; db.yl = 0.0; db.xh = 100.0; db.yh = 20.0;
        for (int i = 0; i < 3; ++i) {
            Node node;
            node.id = i;
            node.name = "warm" + std::to_string(i);
            node.width = 2.0;
            node.height = 2.0;
            node.x = 10.0 + i;
            node.y = 5.0;
            node.fixed = i == 2;
            db.node_by_name.emplace(node.name, node.id);
            db.nodes.push_back(node);
            (node.fixed ? db.fixed_ids : db.movable_ids).push_back(i);
        }
        const Real fixed_x = db.nodes[2].x;
        const std::filesystem::path path =
            std::filesystem::temp_directory_path() / "dpcpp_initial_placement_test.pl";
        {
            std::ofstream out(path);
            out << "UCLA pl 1.0\n"
                << "warm0 20 3 : N\n"
                << "warm1 30 7 : N\n"
                << "warm2 80 10 : N /FIXED\n";
        }
        load_movable_placement(db, path.string());
        require(std::abs(db.nodes[0].x - 21.0) < 1.0e-12 &&
                std::abs(db.nodes[1].y - 8.0) < 1.0e-12,
                "movable warm start failed");
        require(db.nodes[2].x == fixed_x,
                "warm start must not overwrite fixed nodes");
        std::filesystem::remove(path);
    }
    {
        std::vector<Real> values(32);
        for (int i = 0; i < 32; ++i) values[i] = std::sin(0.3 * i) + 0.1 * i;
        const auto expected = values;
        dct2_orthonormal(values, 8, 4);
        idct2_orthonormal(values, 8, 4);
        Real error = 0.0;
        for (int i = 0; i < 32; ++i) error = std::max(error, std::abs(values[i] - expected[i]));
        require(error < 1.0e-10, "DCT round trip failed");
    }
    {
        constexpr int nx = 8;
        constexpr int ny = 4;
        std::vector<Real> values(nx * ny, 0.0);
        values[1] = 1.0;
        inverse_mixed_sine_cosine2(values, nx, ny, true);
        Real error = 0.0;
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const Real expected = std::sqrt(2.0 / nx) *
                    std::sin(3.14159265358979323846 * (x + 0.5) / nx) *
                    std::sqrt(1.0 / ny);
                error = std::max(error, std::abs(values[y * nx + x] - expected));
            }
        }
        require(error < 1.0e-10, "mixed sine/cosine inverse failed");
    }
    {
        Database db;
        db.nodes.resize(3);
        for (int i = 0; i < 3; ++i) {
            db.nodes[i].id = i;
            db.nodes[i].width = db.nodes[i].height = 1.0;
        }
        db.nodes[0].x = 1.0; db.nodes[0].y = 2.0;
        db.nodes[1].x = 6.0; db.nodes[1].y = 4.0;
        db.nodes[2].x = 9.0; db.nodes[2].y = 1.0; db.nodes[2].fixed = true;
        Net net;
        net.pins = {{0, 0.25, -0.1}, {1, -0.4, 0.3}, {2, 0.0, 0.0}};
        db.nets.push_back(net);
        std::vector<Real> gx, gy;
        const Real gamma = 2.0;
        weighted_average_wirelength(db, gamma, 100, &gx, &gy);
        const Real epsilon = 1.0e-5;
        const Real original = db.nodes[0].x;
        db.nodes[0].x = original + epsilon;
        const Real plus = weighted_average_wirelength(db, gamma, 100, nullptr, nullptr);
        db.nodes[0].x = original - epsilon;
        const Real minus = weighted_average_wirelength(db, gamma, 100, nullptr, nullptr);
        db.nodes[0].x = original;
        const Real numeric = (plus - minus) / (2.0 * epsilon);
        require(std::abs(gx[0] - numeric) < 1.0e-6, "WA gradient failed");
        require(gx[2] == 0.0 && gy[2] == 0.0, "fixed-node gradient must be zero");

        std::vector<Real> exact_gx, exact_gy;
        const Real exact_value = exact_hpwl_subgradient(db, 100, &exact_gx, &exact_gy);
        db.nodes[0].x = original + epsilon;
        const Real exact_plus = exact_hpwl(db);
        db.nodes[0].x = original - epsilon;
        const Real exact_minus = exact_hpwl(db);
        db.nodes[0].x = original;
        const Real exact_numeric = (exact_plus - exact_minus) / (2.0 * epsilon);
        require(std::abs(exact_value - exact_hpwl(db)) < 1.0e-12,
                "exact HPWL oracle value mismatch");
        require(std::abs(exact_gx[0] - exact_numeric) < 1.0e-8,
                "exact HPWL pin-offset subgradient failed");
        require(exact_gx[2] == 0.0 && exact_gy[2] == 0.0,
                "exact HPWL fixed-node subgradient must be zero");

        std::vector<Real> group_hpwl;
        std::vector<std::vector<Real>> group_gx, group_gy;
        const Real grouped_value = exact_hpwl_group_subgradients(
            db, 100, 2, group_hpwl, group_gx, group_gy);
        require(std::abs(grouped_value - exact_value) < 1.0e-12 &&
                std::abs(group_hpwl[0] + group_hpwl[1] - exact_value) < 1.0e-12,
                "grouped HPWL value mismatch");
        for (int id = 0; id < 3; ++id) {
            require(std::abs(group_gx[0][id] + group_gx[1][id] - exact_gx[id]) < 1.0e-12 &&
                    std::abs(group_gy[0][id] + group_gy[1][id] - exact_gy[id]) < 1.0e-12,
                    "grouped HPWL subgradient mismatch");
        }
    }
    {
        Database db;
        db.nodes.resize(4);
        for (int i = 0; i < 4; ++i) {
            db.nodes[i].id = i;
            db.nodes[i].x = i < 2 ? 0.0 : 5.0;
            db.nodes[i].y = static_cast<Real>(i);
        }
        Net net;
        net.pins = {{0, 0.0, 0.0}, {1, 0.0, 0.0},
                    {2, 0.0, 0.0}, {3, 0.0, 0.0}};
        db.nets.push_back(net);
        std::vector<Real> gx, gy;
        exact_hpwl_subgradient(db, 100, &gx, &gy);
        require(std::abs(gx[0] + 0.5) < 1.0e-12 &&
                std::abs(gx[1] + 0.5) < 1.0e-12 &&
                std::abs(gx[2] - 0.5) < 1.0e-12 &&
                std::abs(gx[3] - 0.5) < 1.0e-12,
                "exact HPWL tied extrema must split the subgradient");
    }
    {
        Database db;
        db.xl = 0.0; db.yl = 0.0; db.xh = 100.0; db.yh = 20.0;
        db.rows = {{0.0, 10.0, 1.0, 1.0, 0.0, 100},
                   {10.0, 10.0, 1.0, 1.0, 0.0, 100}};
        for (int i = 0; i < 5; ++i) {
            Node node;
            node.id = i;
            node.name = "n" + std::to_string(i);
            node.width = i == 4 ? 20.0 : 15.0;
            node.height = i == 4 ? 20.0 : 10.0;
            node.x = i == 4 ? 50.0 : 48.0 + i;
            node.y = i == 4 ? 10.0 : (i % 2 ? 15.0 : 5.0);
            node.fixed = i == 4;
            db.nodes.push_back(node);
            (node.fixed ? db.fixed_ids : db.movable_ids).push_back(i);
        }
        db.node_pin_weight.assign(db.nodes.size(), 1);
        const LegalizationProxy before_proxy = evaluate_legalization_proxy(db);
        require(std::isfinite(before_proxy.mean_row_distance) &&
                std::isfinite(before_proxy.segment_overflow),
                "legalization proxy must remain finite");
        LegalizeConfig config;
        config.run_dreamplace_detailed = true;
        config.detailed_passes = 1;
        config.detailed_outer_rounds = 2;
        config.cell_insertion_passes = 1;
        config.projected_subgradient_passes = 1;
        config.constrained_bundle_passes = 2;
        config.constrained_bundle_size = 2;
        config.row_relegalization_passes = 1;
        config.independent_set_size = 8;
        config.use_hungarian_matching = true;
        const LegalizeResult result = legalize_and_refine(db, config);
        require(result.legality.legal, "obstacle-aware legalization failed");
        const LegalizationProxy after_proxy = evaluate_legalization_proxy(db);
        require(after_proxy.mean_row_distance < 1.0e-12,
                "legal placement must have zero row-distance proxy");
        for (int id : db.movable_ids) {
            const Node& node = db.nodes[id];
            require(node.x + 0.5 * node.width <= 40.0 + 1.0e-9 ||
                    node.x - 0.5 * node.width >= 60.0 - 1.0e-9,
                    "movable cell overlaps fixed macro");
        }
    }
    std::cout << "core_test passed\n";
    return 0;
}
