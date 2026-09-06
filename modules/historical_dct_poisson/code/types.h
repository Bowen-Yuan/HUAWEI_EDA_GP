#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace dpcpp {

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
    Real area() const { return width * height; }
};

struct Pin {
    int node = -1;
    Real offset_x = 0.0;
    Real offset_y = 0.0;
};

struct Net {
    int id = -1;
    std::string name;
    std::vector<Pin> pins;
    Real weight = 1.0;
};

struct Row {
    Real y = 0.0;
    Real height = 0.0;
    Real site_width = 1.0;
    Real site_spacing = 1.0;
    Real origin = 0.0;
    int num_sites = 0;
    Real xh() const { return origin + num_sites * site_spacing; }
};

struct Database {
    std::vector<Node> nodes;
    std::vector<Net> nets;
    std::vector<Row> rows;
    std::vector<int> movable_ids;
    std::vector<int> fixed_ids;
    std::vector<int> node_pin_weight;
    std::unordered_map<std::string, int> node_by_name;
    Real xl = 0.0;
    Real yl = 0.0;
    Real xh = 0.0;
    Real yh = 0.0;
    std::string benchmark_base;
    std::string raw_pl_path;
};

struct Metrics {
    int iteration = 0;
    Real exact_hpwl = 0.0;
    Real smooth_wirelength = 0.0;
    Real density_energy = 0.0;
    Real overflow = 0.0;
    Real max_density = 0.0;
    Real density_weight = 0.0;
    Real gamma = 0.0;
    Real step = 0.0;
};

}  // namespace dpcpp

