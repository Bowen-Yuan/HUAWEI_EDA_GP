#include "wirelength.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dpcpp {

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

Real wa_axis(const Database& db, const Net& net, Real gamma, bool x_axis,
             std::vector<Real>* gradient) {
    Real minimum = std::numeric_limits<Real>::infinity();
    Real maximum = -std::numeric_limits<Real>::infinity();
    for (const Pin& pin : net.pins) {
        const Node& node = db.nodes[pin.node];
        const Real value = (x_axis ? node.x + pin.offset_x : node.y + pin.offset_y);
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }
    Real bp = 0.0, cp = 0.0, bn = 0.0, cn = 0.0;
    for (const Pin& pin : net.pins) {
        const Node& node = db.nodes[pin.node];
        const Real value = (x_axis ? node.x + pin.offset_x : node.y + pin.offset_y);
        const Real ap = std::exp((value - maximum) / gamma);
        const Real an = std::exp((minimum - value) / gamma);
        bp += ap;
        cp += value * ap;
        bn += an;
        cn += value * an;
    }
    const Real positive_mean = cp / bp;
    const Real negative_mean = cn / bn;
    if (gradient) {
        for (const Pin& pin : net.pins) {
            const Node& node = db.nodes[pin.node];
            if (node.fixed) continue;
            const Real value = (x_axis ? node.x + pin.offset_x : node.y + pin.offset_y);
            const Real ap = std::exp((value - maximum) / gamma);
            const Real an = std::exp((minimum - value) / gamma);
            const Real gp = (ap / bp) * (1.0 + (value - positive_mean) / gamma);
            const Real gn = (an / bn) * (1.0 - (value - negative_mean) / gamma);
            (*gradient)[pin.node] += net.weight * (gp - gn);
        }
    }
    return net.weight * (positive_mean - negative_mean);
}

}  // namespace

Real weighted_average_wirelength(const Database& db, Real gamma,
                                 int degree_limit,
                                 std::vector<Real>* grad_x,
                                 std::vector<Real>* grad_y) {
    gamma = std::max(gamma, 1.0e-9);
    if (grad_x) grad_x->assign(db.nodes.size(), 0.0);
    if (grad_y) grad_y->assign(db.nodes.size(), 0.0);
    Real total = 0.0;
    // Gradient accumulation is intentionally deterministic.  The per-net CPU
    // implementation follows DREAMPlace Algorithm 2 without global atomics.
    for (const Net& net : db.nets) {
        if (degree_limit > 0 && static_cast<int>(net.pins.size()) > degree_limit) continue;
        total += wa_axis(db, net, gamma, true, grad_x);
        total += wa_axis(db, net, gamma, false, grad_y);
    }
    return total;
}

}  // namespace dpcpp
