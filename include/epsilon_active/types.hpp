#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace ea {

using Real = double;

struct Node {
    int id = -1;
    std::string name;
    Real x = 0.0;
    Real y = 0.0;
    Real width = 0.0;
    Real height = 0.0;
    bool fixed = false;
    bool terminal_ni = false;
    std::string orientation = "N";

    Real area() const noexcept { return width * height; }
};

struct Pin {
    int node = -1;
    Real offset_x = 0.0;
    Real offset_y = 0.0;
};

struct Net {
    int id = -1;
    std::string name;
    std::size_t pin_begin = 0;
    std::size_t pin_count = 0;
    Real weight = 1.0;
};

struct Row {
    Real y = 0.0;
    Real height = 0.0;
    Real site_spacing = 1.0;
    Real origin = 0.0;
    int num_sites = 0;

    Real right() const noexcept { return origin + num_sites * site_spacing; }
};

struct Database {
    std::vector<Node> nodes;
    std::vector<Pin> pins;
    std::vector<Net> nets;
    std::vector<Row> rows;
    std::vector<int> movable_ids;
    std::vector<int> fixed_ids;
    std::vector<int> node_pin_count;
    std::vector<std::size_t> node_pin_offsets;
    std::vector<std::size_t> node_pin_indices;
    std::unordered_map<std::string, int> node_by_name;
    Real xl = 0.0;
    Real yl = 0.0;
    Real xh = 0.0;
    Real yh = 0.0;
    Real movable_area = 0.0;
    std::string benchmark_base;
    std::string raw_pl_path;
};

struct DensityMetrics {
    Real energy = 0.0;
    Real overflow = 0.0;
    Real max_density = 0.0;
};

struct IterationMetrics {
    int iteration = 0;
    Real hpwl = 0.0;
    Real overflow = 0.0;
    Real max_density = 0.0;
    Real density_energy = 0.0;
    Real lambda_base = 0.0;
    Real lambda_control = 1.0;
    Real lambda_effective = 0.0;
    Real hpwl_epsilon = 0.0;
    Real density_epsilon = 0.0;
    Real epsilon_scale = 1.0;
    Real learning_rate = 0.0;
};

}  // namespace ea

