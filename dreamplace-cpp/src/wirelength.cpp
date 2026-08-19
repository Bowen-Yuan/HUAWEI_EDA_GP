#include "wirelength.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dpcpp {

namespace {

struct ActiveWorkspace {
    const Database* database = nullptr;
    std::vector<std::size_t> pin_offsets;
    std::vector<std::size_t> node_offsets;
    std::vector<std::size_t> node_pin_indices;
    std::vector<int> pin_nodes;
    std::vector<Real> pin_offset_x;
    std::vector<Real> pin_offset_y;
    std::vector<Real> net_weights;
    std::vector<Real> net_values;
    std::vector<Real> pin_gx;
    std::vector<Real> pin_gy;
};

void prepare_active_topology(const Database& db, ActiveWorkspace& workspace) {
    if (workspace.database == &db &&
        workspace.pin_offsets.size() == db.nets.size() + 1 &&
        workspace.node_offsets.size() == db.nodes.size() + 1) return;
    workspace.database = &db;
    workspace.pin_offsets.resize(db.nets.size() + 1);
    workspace.pin_offsets[0] = 0;
    for (std::size_t n = 0; n < db.nets.size(); ++n) {
        workspace.pin_offsets[n + 1] =
            workspace.pin_offsets[n] + db.nets[n].pins.size();
    }
    workspace.node_offsets.assign(db.nodes.size() + 1, 0);
    for (const Net& net : db.nets) {
        for (const Pin& pin : net.pins)
            ++workspace.node_offsets[pin.node + 1];
    }
    for (std::size_t i = 1; i < workspace.node_offsets.size(); ++i)
        workspace.node_offsets[i] += workspace.node_offsets[i - 1];
    workspace.node_pin_indices.resize(workspace.pin_offsets.back());
    workspace.pin_nodes.resize(workspace.pin_offsets.back());
    workspace.pin_offset_x.resize(workspace.pin_offsets.back());
    workspace.pin_offset_y.resize(workspace.pin_offsets.back());
    workspace.net_weights.resize(db.nets.size());
    std::vector<std::size_t> cursor = workspace.node_offsets;
    for (std::size_t n = 0; n < db.nets.size(); ++n) {
        const std::size_t offset = workspace.pin_offsets[n];
        workspace.net_weights[n] = db.nets[n].weight;
        for (std::size_t p = 0; p < db.nets[n].pins.size(); ++p) {
            const Pin& source = db.nets[n].pins[p];
            const int node = source.node;
            workspace.pin_nodes[offset + p] = node;
            workspace.pin_offset_x[offset + p] = source.offset_x;
            workspace.pin_offset_y[offset + p] = source.offset_y;
            workspace.node_pin_indices[cursor[node]++] = offset + p;
        }
    }
}

inline Real active_weight(Real value, Real power) {
    if (value <= 0.0) return 0.0;
    if (power == 1.0) return value;
    const Real square = value * value;
    if (power == 2.0) return square;
    if (power == 4.0) return square * square;
    return std::pow(value, power);
}

}  // namespace

Real exact_hpwl(const Database& db, int degree_limit) {
    Real total = 0.0;
    #pragma omp parallel for reduction(+:total) schedule(dynamic, 512)
    for (int n = 0; n < static_cast<int>(db.nets.size()); ++n) {
        const Net& net = db.nets[n];
        if (degree_limit > 0 && static_cast<int>(net.pins.size()) > degree_limit) continue;
        Real min_x = std::numeric_limits<Real>::infinity();
        Real max_x = -std::numeric_limits<Real>::infinity();
        Real min_y = std::numeric_limits<Real>::infinity();
        Real max_y = -std::numeric_limits<Real>::infinity();
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            min_x = std::min(min_x, x);
            max_x = std::max(max_x, x);
            min_y = std::min(min_y, y);
            max_y = std::max(max_y, y);
        }
        total += net.weight * ((max_x - min_x) + (max_y - min_y));
    }
    return total;
}

Real exact_hpwl_subgradient(const Database& db, int gradient_degree_limit,
                            std::vector<Real>* grad_x,
                            std::vector<Real>* grad_y) {
    if (grad_x) grad_x->assign(db.nodes.size(), 0.0);
    if (grad_y) grad_y->assign(db.nodes.size(), 0.0);
    Real total = 0.0;
    for (const Net& net : db.nets) {
        if (net.pins.size() < 2) continue;
        Real min_x = std::numeric_limits<Real>::infinity();
        Real max_x = -std::numeric_limits<Real>::infinity();
        Real min_y = std::numeric_limits<Real>::infinity();
        Real max_y = -std::numeric_limits<Real>::infinity();
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            min_x = std::min(min_x, x); max_x = std::max(max_x, x);
            min_y = std::min(min_y, y); max_y = std::max(max_y, y);
        }
        total += net.weight * ((max_x - min_x) + (max_y - min_y));
        if ((!grad_x && !grad_y) ||
            (gradient_degree_limit > 0 &&
             static_cast<int>(net.pins.size()) > gradient_degree_limit)) continue;

        int count_min_x = 0, count_max_x = 0;
        int count_min_y = 0, count_max_y = 0;
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            if (x == min_x) ++count_min_x;
            if (x == max_x) ++count_max_x;
            if (y == min_y) ++count_min_y;
            if (y == max_y) ++count_max_y;
        }
        const Real min_x_share = net.weight / std::max(1, count_min_x);
        const Real max_x_share = net.weight / std::max(1, count_max_x);
        const Real min_y_share = net.weight / std::max(1, count_min_y);
        const Real max_y_share = net.weight / std::max(1, count_max_y);
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            if (node.fixed) continue;
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            if (grad_x) {
                if (x == max_x) (*grad_x)[pin.node] += max_x_share;
                if (x == min_x) (*grad_x)[pin.node] -= min_x_share;
            }
            if (grad_y) {
                if (y == max_y) (*grad_y)[pin.node] += max_y_share;
                if (y == min_y) (*grad_y)[pin.node] -= min_y_share;
            }
        }
    }
    return total;
}

Real exact_hpwl_active_set_direction(const Database& db,
                                     int gradient_degree_limit,
                                     Real radius,
                                     Real power,
                                     std::vector<Real>* grad_x,
                                     std::vector<Real>* grad_y,
                                     Real span_cap,
                                     Real small_span_min_radius,
                                     Real small_span_threshold,
                                     Real span_cap_blend) {
    if (radius <= 0.0)
        return exact_hpwl_subgradient(
            db, gradient_degree_limit, grad_x, grad_y);
    if (grad_x) grad_x->assign(db.nodes.size(), 0.0);
    if (grad_y) grad_y->assign(db.nodes.size(), 0.0);
    static ActiveWorkspace workspace;
    prepare_active_topology(db, workspace);
    workspace.net_values.resize(db.nets.size());
    workspace.pin_gx.resize(workspace.pin_offsets.back());
    workspace.pin_gy.resize(workspace.pin_offsets.back());
    #pragma omp parallel for schedule(static, 256)
    for (int n = 0; n < static_cast<int>(db.nets.size()); ++n) {
        const std::size_t offset = workspace.pin_offsets[n];
        const std::size_t end = workspace.pin_offsets[n + 1];
        const std::size_t pin_count = end - offset;
        if (pin_count < 2) {
            workspace.net_values[n] = 0.0;
            continue;
        }
        Real min_x = std::numeric_limits<Real>::infinity();
        Real max_x = -std::numeric_limits<Real>::infinity();
        Real min_y = std::numeric_limits<Real>::infinity();
        Real max_y = -std::numeric_limits<Real>::infinity();
        for (std::size_t pin_index = offset; pin_index < end; ++pin_index) {
            const int node_id = workspace.pin_nodes[pin_index];
            const Node& node = db.nodes[node_id];
            const Real x = node.x + workspace.pin_offset_x[pin_index];
            const Real y = node.y + workspace.pin_offset_y[pin_index];
            min_x = std::min(min_x, x); max_x = std::max(max_x, x);
            min_y = std::min(min_y, y); max_y = std::max(max_y, y);
        }
        workspace.net_values[n] =
            workspace.net_weights[n] * ((max_x - min_x) + (max_y - min_y));
        if ((!grad_x && !grad_y) ||
            (gradient_degree_limit > 0 && static_cast<int>(pin_count) >
                                              gradient_degree_limit)) {
            std::fill(workspace.pin_gx.begin() + offset,
                      workspace.pin_gx.begin() + end, 0.0);
            std::fill(workspace.pin_gy.begin() + offset,
                      workspace.pin_gy.begin() + end, 0.0);
            continue;
        }

        const auto per_axis_radius = [&](Real span) {
            if (span_cap <= 0.0) return radius;
            Real local_radius = span_cap * span;
            if (small_span_threshold > 0.0 && span < small_span_threshold) {
                local_radius = std::max(local_radius, small_span_min_radius);
            }
            const Real relative_radius = std::min(
                radius, std::max<Real>(1.0e-12, local_radius));
            return (1.0 - span_cap_blend) * radius +
                   span_cap_blend * relative_radius;
        };
        const Real radius_x = per_axis_radius(max_x - min_x);
        const Real radius_y = per_axis_radius(max_y - min_y);
        Real min_x_weight = 0.0, max_x_weight = 0.0;
        Real min_y_weight = 0.0, max_y_weight = 0.0;
        for (std::size_t pin_index = offset; pin_index < end; ++pin_index) {
            const int node_id = workspace.pin_nodes[pin_index];
            const Node& node = db.nodes[node_id];
            const Real x = node.x + workspace.pin_offset_x[pin_index];
            const Real y = node.y + workspace.pin_offset_y[pin_index];
            const Real wx_min = std::max<Real>(0.0, 1.0 - (x - min_x) / radius_x);
            const Real wx_max = std::max<Real>(0.0, 1.0 - (max_x - x) / radius_x);
            const Real wy_min = std::max<Real>(0.0, 1.0 - (y - min_y) / radius_y);
            const Real wy_max = std::max<Real>(0.0, 1.0 - (max_y - y) / radius_y);
            min_x_weight += active_weight(wx_min, power);
            max_x_weight += active_weight(wx_max, power);
            min_y_weight += active_weight(wy_min, power);
            max_y_weight += active_weight(wy_max, power);
        }
        for (std::size_t pin_index = offset; pin_index < end; ++pin_index) {
            const int node_id = workspace.pin_nodes[pin_index];
            const Node& node = db.nodes[node_id];
            if (node.fixed) {
                workspace.pin_gx[pin_index] = 0.0;
                workspace.pin_gy[pin_index] = 0.0;
                continue;
            }
            const Real x = node.x + workspace.pin_offset_x[pin_index];
            const Real y = node.y + workspace.pin_offset_y[pin_index];
            if (grad_x || grad_y) {
                const Real wx_min = std::max<Real>(0.0, 1.0 - (x - min_x) / radius_x);
                const Real wx_max = std::max<Real>(0.0, 1.0 - (max_x - x) / radius_x);
                workspace.pin_gx[pin_index] = workspace.net_weights[n] *
                    (active_weight(wx_max, power) /
                         std::max<Real>(max_x_weight, 1.0e-30) -
                     active_weight(wx_min, power) /
                         std::max<Real>(min_x_weight, 1.0e-30));
                const Real wy_min = std::max<Real>(0.0, 1.0 - (y - min_y) / radius_y);
                const Real wy_max = std::max<Real>(0.0, 1.0 - (max_y - y) / radius_y);
                workspace.pin_gy[pin_index] = workspace.net_weights[n] *
                    (active_weight(wy_max, power) /
                         std::max<Real>(max_y_weight, 1.0e-30) -
                     active_weight(wy_min, power) /
                         std::max<Real>(min_y_weight, 1.0e-30));
            }
        }
    }
    Real total = 0.0;
    #pragma omp parallel for reduction(+:total) schedule(static, 256)
    for (int n = 0; n < static_cast<int>(workspace.net_values.size()); ++n)
        total += workspace.net_values[n];
    if (!grad_x && !grad_y) return total;
    #pragma omp parallel for schedule(static)
    for (int node = 0; node < static_cast<int>(db.nodes.size()); ++node) {
        if (db.nodes[node].fixed) continue;
        Real sum_x = 0.0;
        Real sum_y = 0.0;
        for (std::size_t i = workspace.node_offsets[node];
             i < workspace.node_offsets[node + 1]; ++i) {
            const std::size_t pin = workspace.node_pin_indices[i];
            sum_x += workspace.pin_gx[pin];
            sum_y += workspace.pin_gy[pin];
        }
        if (grad_x) (*grad_x)[node] = sum_x;
        if (grad_y) (*grad_y)[node] = sum_y;
    }
    return total;
}

Real exact_hpwl_primal_dual_direction(const Database& db,
                                      int gradient_degree_limit,
                                      Real dual_step,
                                      std::vector<Real>* grad_x,
                                      std::vector<Real>* grad_y) {
    if (dual_step <= 0.0)
        return exact_hpwl_subgradient(
            db, gradient_degree_limit, grad_x, grad_y);
    if (grad_x) grad_x->assign(db.nodes.size(), 0.0);
    if (grad_y) grad_y->assign(db.nodes.size(), 0.0);

    struct DualWorkspace {
        const Database* database = nullptr;
        ActiveWorkspace topology;
        std::vector<Real> positive_x;
        std::vector<Real> negative_x;
        std::vector<Real> positive_y;
        std::vector<Real> negative_y;
    };
    static DualWorkspace dual;
    prepare_active_topology(db, dual.topology);
    const std::size_t total_pins = dual.topology.pin_offsets.back();
    if (dual.database != &db || dual.positive_x.size() != total_pins) {
        dual.database = &db;
        dual.positive_x.assign(total_pins, 0.0);
        dual.negative_x.assign(total_pins, 0.0);
        dual.positive_y.assign(total_pins, 0.0);
        dual.negative_y.assign(total_pins, 0.0);
        for (std::size_t n = 0; n < db.nets.size(); ++n) {
            const std::size_t begin = dual.topology.pin_offsets[n];
            const std::size_t end = dual.topology.pin_offsets[n + 1];
            const Real uniform = 1.0 /
                std::max<Real>(1.0, static_cast<Real>(end - begin));
            std::fill(dual.positive_x.begin() + begin,
                      dual.positive_x.begin() + end, uniform);
            std::fill(dual.negative_x.begin() + begin,
                      dual.negative_x.begin() + end, uniform);
            std::fill(dual.positive_y.begin() + begin,
                      dual.positive_y.begin() + end, uniform);
            std::fill(dual.negative_y.begin() + begin,
                      dual.negative_y.begin() + end, uniform);
        }
    }

    ActiveWorkspace& workspace = dual.topology;
    workspace.net_values.resize(db.nets.size());
    workspace.pin_gx.resize(total_pins);
    workspace.pin_gy.resize(total_pins);
    const Real inv_width = 1.0 / std::max<Real>(db.xh - db.xl, 1.0);
    const Real inv_height = 1.0 / std::max<Real>(db.yh - db.yl, 1.0);

    #pragma omp parallel for schedule(static, 256)
    for (int n = 0; n < static_cast<int>(db.nets.size()); ++n) {
        const Net& net = db.nets[n];
        const std::size_t begin = workspace.pin_offsets[n];
        const std::size_t degree = net.pins.size();
        Real min_x = std::numeric_limits<Real>::infinity();
        Real max_x = -std::numeric_limits<Real>::infinity();
        Real min_y = std::numeric_limits<Real>::infinity();
        Real max_y = -std::numeric_limits<Real>::infinity();
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            min_x = std::min(min_x, x); max_x = std::max(max_x, x);
            min_y = std::min(min_y, y); max_y = std::max(max_y, y);
        }
        workspace.net_values[n] = degree >= 2
            ? net.weight * ((max_x - min_x) + (max_y - min_y))
            : 0.0;
        if (degree < 2 || (gradient_degree_limit > 0 &&
                           static_cast<int>(degree) > gradient_degree_limit)) {
            std::fill(workspace.pin_gx.begin() + begin,
                      workspace.pin_gx.begin() + begin + degree, 0.0);
            std::fill(workspace.pin_gy.begin() + begin,
                      workspace.pin_gy.begin() + begin + degree, 0.0);
            continue;
        }

        std::vector<Real> scratch(degree);
        auto project_simplex = [&](Real* values) {
            std::copy(values, values + degree, scratch.begin());
            std::sort(scratch.begin(), scratch.end(), std::greater<Real>());
            Real prefix = 0.0;
            Real theta = 0.0;
            for (std::size_t i = 0; i < degree; ++i) {
                prefix += scratch[i];
                const Real candidate = (prefix - 1.0) /
                    static_cast<Real>(i + 1);
                if (scratch[i] > candidate) theta = candidate;
            }
            for (std::size_t i = 0; i < degree; ++i)
                values[i] = std::max<Real>(0.0, values[i] - theta);
        };

        for (std::size_t p = 0; p < degree; ++p) {
            const Pin& pin = net.pins[p];
            const Node& node = db.nodes[pin.node];
            const Real x = (node.x + pin.offset_x - db.xl) * inv_width;
            const Real y = (node.y + pin.offset_y - db.yl) * inv_height;
            dual.positive_x[begin + p] += dual_step * x;
            dual.negative_x[begin + p] -= dual_step * x;
            dual.positive_y[begin + p] += dual_step * y;
            dual.negative_y[begin + p] -= dual_step * y;
        }
        project_simplex(dual.positive_x.data() + begin);
        project_simplex(dual.negative_x.data() + begin);
        project_simplex(dual.positive_y.data() + begin);
        project_simplex(dual.negative_y.data() + begin);
        for (std::size_t p = 0; p < degree; ++p) {
            const Pin& pin = net.pins[p];
            if (db.nodes[pin.node].fixed) {
                workspace.pin_gx[begin + p] = 0.0;
                workspace.pin_gy[begin + p] = 0.0;
            } else {
                workspace.pin_gx[begin + p] = net.weight *
                    (dual.positive_x[begin + p] - dual.negative_x[begin + p]);
                workspace.pin_gy[begin + p] = net.weight *
                    (dual.positive_y[begin + p] - dual.negative_y[begin + p]);
            }
        }
    }

    Real total = 0.0;
    for (Real value : workspace.net_values) total += value;
    #pragma omp parallel for schedule(static)
    for (int node = 0; node < static_cast<int>(db.nodes.size()); ++node) {
        if (db.nodes[node].fixed) continue;
        Real sum_x = 0.0;
        Real sum_y = 0.0;
        for (std::size_t i = workspace.node_offsets[node];
             i < workspace.node_offsets[node + 1]; ++i) {
            const std::size_t pin = workspace.node_pin_indices[i];
            sum_x += workspace.pin_gx[pin];
            sum_y += workspace.pin_gy[pin];
        }
        if (grad_x) (*grad_x)[node] = sum_x;
        if (grad_y) (*grad_y)[node] = sum_y;
    }
    return total;
}

Real exact_hpwl_group_subgradients(
    const Database& db, int gradient_degree_limit, int group_count,
    std::vector<Real>& group_hpwl,
    std::vector<std::vector<Real>>& group_grad_x,
    std::vector<std::vector<Real>>& group_grad_y) {
    if (group_count <= 0) throw std::runtime_error("HPWL group count must be positive");
    group_hpwl.assign(group_count, 0.0);
    group_grad_x.assign(group_count, std::vector<Real>(db.nodes.size(), 0.0));
    group_grad_y.assign(group_count, std::vector<Real>(db.nodes.size(), 0.0));
    Real total = 0.0;
    for (int net_index = 0; net_index < static_cast<int>(db.nets.size()); ++net_index) {
        const Net& net = db.nets[net_index];
        if (net.pins.size() < 2) continue;
        const int group = net_index % group_count;
        Real min_x = std::numeric_limits<Real>::infinity();
        Real max_x = -std::numeric_limits<Real>::infinity();
        Real min_y = std::numeric_limits<Real>::infinity();
        Real max_y = -std::numeric_limits<Real>::infinity();
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            min_x = std::min(min_x, x); max_x = std::max(max_x, x);
            min_y = std::min(min_y, y); max_y = std::max(max_y, y);
        }
        const Real value = net.weight * ((max_x - min_x) + (max_y - min_y));
        total += value;
        group_hpwl[group] += value;
        if (gradient_degree_limit > 0 &&
            static_cast<int>(net.pins.size()) > gradient_degree_limit) continue;

        int count_min_x = 0, count_max_x = 0;
        int count_min_y = 0, count_max_y = 0;
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            if (x == min_x) ++count_min_x;
            if (x == max_x) ++count_max_x;
            if (y == min_y) ++count_min_y;
            if (y == max_y) ++count_max_y;
        }
        const Real min_x_share = net.weight / std::max(1, count_min_x);
        const Real max_x_share = net.weight / std::max(1, count_max_x);
        const Real min_y_share = net.weight / std::max(1, count_min_y);
        const Real max_y_share = net.weight / std::max(1, count_max_y);
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            if (node.fixed) continue;
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            if (x == max_x) group_grad_x[group][pin.node] += max_x_share;
            if (x == min_x) group_grad_x[group][pin.node] -= min_x_share;
            if (y == max_y) group_grad_y[group][pin.node] += max_y_share;
            if (y == min_y) group_grad_y[group][pin.node] -= min_y_share;
        }
    }
    return total;
}

namespace {

struct WaWorkspace {
    const Database* database = nullptr;
    std::vector<std::size_t> pin_offsets;
    std::vector<std::size_t> node_offsets;
    std::vector<std::size_t> node_pin_indices;
    std::vector<Real> value_x;
    std::vector<Real> value_y;
    std::vector<Real> pin_positive_x;
    std::vector<Real> pin_negative_x;
    std::vector<Real> pin_positive_y;
    std::vector<Real> pin_negative_y;
};

void prepare_wa_topology(const Database& db, WaWorkspace& workspace) {
    if (workspace.database == &db &&
        workspace.pin_offsets.size() == db.nets.size() + 1 &&
        workspace.node_offsets.size() == db.nodes.size() + 1) return;
    workspace.database = &db;
    workspace.pin_offsets.resize(db.nets.size() + 1);
    workspace.pin_offsets[0] = 0;
    for (std::size_t n = 0; n < db.nets.size(); ++n)
        workspace.pin_offsets[n + 1] =
            workspace.pin_offsets[n] + db.nets[n].pins.size();

    workspace.node_offsets.assign(db.nodes.size() + 1, 0);
    for (const Net& net : db.nets) {
        for (const Pin& pin : net.pins)
            ++workspace.node_offsets[pin.node + 1];
    }
    for (std::size_t i = 1; i < workspace.node_offsets.size(); ++i)
        workspace.node_offsets[i] += workspace.node_offsets[i - 1];
    workspace.node_pin_indices.resize(workspace.pin_offsets.back());
    std::vector<std::size_t> cursor = workspace.node_offsets;
    for (std::size_t n = 0; n < db.nets.size(); ++n) {
        const std::size_t offset = workspace.pin_offsets[n];
        for (std::size_t p = 0; p < db.nets[n].pins.size(); ++p) {
            const int node = db.nets[n].pins[p].node;
            workspace.node_pin_indices[cursor[node]++] = offset + p;
        }
    }
}

void wa_net(const Database& db, const Net& net, Real gamma,
            Real& result_x, Real& result_y,
            Real* pin_positive_x, Real* pin_negative_x,
            Real* pin_positive_y, Real* pin_negative_y) {
    Real min_x = std::numeric_limits<Real>::infinity();
    Real max_x = -std::numeric_limits<Real>::infinity();
    Real min_y = std::numeric_limits<Real>::infinity();
    Real max_y = -std::numeric_limits<Real>::infinity();
    for (std::size_t p = 0; p < net.pins.size(); ++p) {
        const Pin& pin = net.pins[p];
        const Node& node = db.nodes[pin.node];
        const Real x = node.x + pin.offset_x;
        const Real y = node.y + pin.offset_y;
        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
    }
    Real bpx = 0.0, cpx = 0.0, bnx = 0.0, cnx = 0.0;
    Real bpy = 0.0, cpy = 0.0, bny = 0.0, cny = 0.0;
    for (std::size_t p = 0; p < net.pins.size(); ++p) {
        const Pin& pin = net.pins[p];
        const Node& node = db.nodes[pin.node];
        const Real x = node.x + pin.offset_x;
        const Real y = node.y + pin.offset_y;
        const Real apx = std::exp((x - max_x) / gamma);
        const Real anx = std::exp((min_x - x) / gamma);
        const Real apy = std::exp((y - max_y) / gamma);
        const Real any = std::exp((min_y - y) / gamma);
        if (pin_positive_x) {
            pin_positive_x[p] = apx;
            pin_negative_x[p] = anx;
            pin_positive_y[p] = apy;
            pin_negative_y[p] = any;
        }
        bpx += apx;
        cpx += x * apx;
        bnx += anx;
        cnx += x * anx;
        bpy += apy;
        cpy += y * apy;
        bny += any;
        cny += y * any;
    }
    const Real positive_x = cpx / bpx;
    const Real negative_x = cnx / bnx;
    const Real positive_y = cpy / bpy;
    const Real negative_y = cny / bny;
    result_x = net.weight * (positive_x - negative_x);
    result_y = net.weight * (positive_y - negative_y);
    if (pin_positive_x) {
        for (std::size_t p = 0; p < net.pins.size(); ++p) {
            const Pin& pin = net.pins[p];
            const Node& node = db.nodes[pin.node];
            if (node.fixed) {
                pin_positive_x[p] = 0.0;
                pin_positive_y[p] = 0.0;
                continue;
            }
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            const Real apx = pin_positive_x[p];
            const Real anx = pin_negative_x[p];
            const Real apy = pin_positive_y[p];
            const Real any = pin_negative_y[p];
            const Real gpx = (apx / bpx) * (1.0 + (x - positive_x) / gamma);
            const Real gnx = (anx / bnx) * (1.0 - (x - negative_x) / gamma);
            const Real gpy = (apy / bpy) * (1.0 + (y - positive_y) / gamma);
            const Real gny = (any / bny) * (1.0 - (y - negative_y) / gamma);
            pin_positive_x[p] = net.weight * (gpx - gnx);
            pin_positive_y[p] = net.weight * (gpy - gny);
        }
    }
}

}  // namespace

Real weighted_average_wirelength(const Database& db, Real gamma,
                                 int degree_limit,
                                 std::vector<Real>* grad_x,
                                 std::vector<Real>* grad_y) {
    gamma = std::max(gamma, 1.0e-9);
    if (grad_x) grad_x->assign(db.nodes.size(), 0.0);
    if (grad_y) grad_y->assign(db.nodes.size(), 0.0);
    // Placement invokes this kernel serially; the OpenMP region below only
    // parallelizes nets within one invocation. Reusing the buffers avoids
    // allocating tens of megabytes twice per optimizer iteration.
    static WaWorkspace workspace;
    prepare_wa_topology(db, workspace);
    workspace.value_x.resize(db.nets.size());
    workspace.value_y.resize(db.nets.size());
    const bool need_gradient = grad_x || grad_y;
    if (need_gradient) {
        workspace.pin_positive_x.resize(workspace.pin_offsets.back());
        workspace.pin_negative_x.resize(workspace.pin_offsets.back());
        workspace.pin_positive_y.resize(workspace.pin_offsets.back());
        workspace.pin_negative_y.resize(workspace.pin_offsets.back());
    }

    #pragma omp parallel for schedule(static, 256)
    for (int n = 0; n < static_cast<int>(db.nets.size()); ++n) {
        const Net& net = db.nets[n];
        if (degree_limit > 0 && static_cast<int>(net.pins.size()) > degree_limit) {
            workspace.value_x[n] = 0.0;
            workspace.value_y[n] = 0.0;
            if (need_gradient) {
                const std::size_t begin = workspace.pin_offsets[n];
                const std::size_t end = workspace.pin_offsets[n + 1];
                std::fill(workspace.pin_positive_x.begin() + begin,
                          workspace.pin_positive_x.begin() + end, 0.0);
                std::fill(workspace.pin_positive_y.begin() + begin,
                          workspace.pin_positive_y.begin() + end, 0.0);
            }
            continue;
        }
        const std::size_t offset = workspace.pin_offsets[n];
        wa_net(db, net, gamma, workspace.value_x[n], workspace.value_y[n],
               need_gradient ? workspace.pin_positive_x.data() + offset : nullptr,
               need_gradient ? workspace.pin_negative_x.data() + offset : nullptr,
               need_gradient ? workspace.pin_positive_y.data() + offset : nullptr,
               need_gradient ? workspace.pin_negative_y.data() + offset : nullptr);
    }

    Real total = 0.0;
    for (std::size_t n = 0; n < db.nets.size(); ++n) {
        const Net& net = db.nets[n];
        if (degree_limit > 0 && static_cast<int>(net.pins.size()) > degree_limit)
            continue;
        total += workspace.value_x[n];
        total += workspace.value_y[n];
    }
    if (!need_gradient) return total;
    #pragma omp parallel for schedule(static)
    for (int node = 0; node < static_cast<int>(db.nodes.size()); ++node) {
        if (db.nodes[node].fixed) continue;
        Real sum_x = 0.0;
        Real sum_y = 0.0;
        for (std::size_t i = workspace.node_offsets[node];
             i < workspace.node_offsets[node + 1]; ++i) {
            const std::size_t pin = workspace.node_pin_indices[i];
            sum_x += workspace.pin_positive_x[pin];
            sum_y += workspace.pin_positive_y[pin];
        }
        if (grad_x) (*grad_x)[node] = sum_x;
        if (grad_y) (*grad_y)[node] = sum_y;
    }
    return total;
}

}  // namespace dpcpp
