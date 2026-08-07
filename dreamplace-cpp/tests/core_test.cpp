#include "bookshelf.h"
#include "legalizer.h"
#include "spectral.h"
#include "wirelength.h"

#include <algorithm>
#include <cmath>
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
        LegalizeConfig config;
        config.run_dreamplace_detailed = true;
        config.detailed_passes = 1;
        const LegalizeResult result = legalize_and_refine(db, config);
        require(result.legality.legal, "obstacle-aware legalization failed");
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
