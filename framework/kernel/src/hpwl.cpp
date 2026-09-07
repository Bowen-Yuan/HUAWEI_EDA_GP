#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ea {
namespace {

Real active_weight(Real value, Real power) {
    if (value <= 0.0) return 0.0;
    if (power == 1.0) return value;
    if (power == 2.0) return value * value;
    if (power == 4.0) {
        const Real square = value * value;
        return square * square;
    }
    return std::pow(value, power);
}

}  // namespace

ExactHpwl::ExactHpwl(const Database& db)
    : db_(db), pin_grad_x_(db.pins.size()), pin_grad_y_(db.pins.size()) {}

Real ExactHpwl::evaluate(Real epsilon, Real active_power, int degree_limit,
                         std::vector<Real>* grad_x,
                         std::vector<Real>* grad_y) {
    if (epsilon < 0.0 || active_power <= 0.0) {
        throw std::invalid_argument("invalid HPWL active-set parameters");
    }
    const bool need_gradient = grad_x != nullptr || grad_y != nullptr;
    Real total = 0.0;
    #pragma omp parallel for reduction(+:total) schedule(dynamic, 256)
    for (int net_index = 0; net_index < static_cast<int>(db_.nets.size()); ++net_index) {
        const Net& net = db_.nets[net_index];
        const std::size_t begin = net.pin_begin;
        const std::size_t end = begin + net.pin_count;
        Real min_x = std::numeric_limits<Real>::infinity();
        Real max_x = -std::numeric_limits<Real>::infinity();
        Real min_y = std::numeric_limits<Real>::infinity();
        Real max_y = -std::numeric_limits<Real>::infinity();
        for (std::size_t p = begin; p < end; ++p) {
            const Pin& pin = db_.pins[p];
            const Node& node = db_.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            min_x = std::min(min_x, x);
            max_x = std::max(max_x, x);
            min_y = std::min(min_y, y);
            max_y = std::max(max_y, y);
        }
        total += net.weight * ((max_x - min_x) + (max_y - min_y));
        if (!need_gradient || (degree_limit > 0 &&
                               static_cast<int>(net.pin_count) > degree_limit)) {
            if (need_gradient) {
                std::fill(pin_grad_x_.begin() + begin, pin_grad_x_.begin() + end, 0.0);
                std::fill(pin_grad_y_.begin() + begin, pin_grad_y_.begin() + end, 0.0);
            }
            continue;
        }

        Real min_x_sum = 0.0;
        Real max_x_sum = 0.0;
        Real min_y_sum = 0.0;
        Real max_y_sum = 0.0;
        for (std::size_t p = begin; p < end; ++p) {
            const Pin& pin = db_.pins[p];
            const Node& node = db_.nodes[pin.node];
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            const Real wx_min = epsilon > 0.0
                ? active_weight(std::max<Real>(0.0, 1.0 - (x - min_x) / epsilon), active_power)
                : (x == min_x ? 1.0 : 0.0);
            const Real wx_max = epsilon > 0.0
                ? active_weight(std::max<Real>(0.0, 1.0 - (max_x - x) / epsilon), active_power)
                : (x == max_x ? 1.0 : 0.0);
            const Real wy_min = epsilon > 0.0
                ? active_weight(std::max<Real>(0.0, 1.0 - (y - min_y) / epsilon), active_power)
                : (y == min_y ? 1.0 : 0.0);
            const Real wy_max = epsilon > 0.0
                ? active_weight(std::max<Real>(0.0, 1.0 - (max_y - y) / epsilon), active_power)
                : (y == max_y ? 1.0 : 0.0);
            min_x_sum += wx_min;
            max_x_sum += wx_max;
            min_y_sum += wy_min;
            max_y_sum += wy_max;
        }
        for (std::size_t p = begin; p < end; ++p) {
            const Pin& pin = db_.pins[p];
            const Node& node = db_.nodes[pin.node];
            if (node.fixed) {
                pin_grad_x_[p] = 0.0;
                pin_grad_y_[p] = 0.0;
                continue;
            }
            const Real x = node.x + pin.offset_x;
            const Real y = node.y + pin.offset_y;
            const Real wx_min = epsilon > 0.0
                ? active_weight(std::max<Real>(0.0, 1.0 - (x - min_x) / epsilon), active_power)
                : (x == min_x ? 1.0 : 0.0);
            const Real wx_max = epsilon > 0.0
                ? active_weight(std::max<Real>(0.0, 1.0 - (max_x - x) / epsilon), active_power)
                : (x == max_x ? 1.0 : 0.0);
            const Real wy_min = epsilon > 0.0
                ? active_weight(std::max<Real>(0.0, 1.0 - (y - min_y) / epsilon), active_power)
                : (y == min_y ? 1.0 : 0.0);
            const Real wy_max = epsilon > 0.0
                ? active_weight(std::max<Real>(0.0, 1.0 - (max_y - y) / epsilon), active_power)
                : (y == max_y ? 1.0 : 0.0);
            pin_grad_x_[p] = net.weight *
                (wx_max / std::max<Real>(max_x_sum, 1.0e-30) -
                 wx_min / std::max<Real>(min_x_sum, 1.0e-30));
            pin_grad_y_[p] = net.weight *
                (wy_max / std::max<Real>(max_y_sum, 1.0e-30) -
                 wy_min / std::max<Real>(min_y_sum, 1.0e-30));
        }
    }

    if (grad_x) grad_x->assign(db_.nodes.size(), 0.0);
    if (grad_y) grad_y->assign(db_.nodes.size(), 0.0);
    if (!need_gradient) return total;
    #pragma omp parallel for schedule(static)
    for (int node_id = 0; node_id < static_cast<int>(db_.nodes.size()); ++node_id) {
        Real sum_x = 0.0;
        Real sum_y = 0.0;
        for (std::size_t i = db_.node_pin_offsets[node_id];
             i < db_.node_pin_offsets[node_id + 1]; ++i) {
            const std::size_t pin = db_.node_pin_indices[i];
            sum_x += pin_grad_x_[pin];
            sum_y += pin_grad_y_[pin];
        }
        if (grad_x) (*grad_x)[node_id] = sum_x;
        if (grad_y) (*grad_y)[node_id] = sum_y;
    }
    return total;
}

}  // namespace ea

