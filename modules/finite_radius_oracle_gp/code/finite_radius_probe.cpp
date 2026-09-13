#include "finite_radius_probe.hpp"

#include <algorithm>
#include <stdexcept>

namespace nsgp::modules {

ea::AuxiliaryQueryStats finite_radius_probe(
    const ea::Database& db, ea::ExactOverlapDensity& density,
    const ea::AuxiliaryContext& context,
    int probe_count, ea::Real radius_bin_scale, ea::Real descent_tolerance,
    std::vector<ea::Real>& aux_x, std::vector<ea::Real>& aux_y) {
    if (probe_count < 0 || radius_bin_scale <= 0.0 || descent_tolerance < 0.0) {
        throw std::invalid_argument(
            "finite_radius_probe requires probe_count >= 0, "
            "radius_bin_scale > 0, descent_tolerance >= 0");
    }
    aux_x.assign(db.nodes.size(), 0.0);
    aux_y.assign(db.nodes.size(), 0.0);
    ea::AuxiliaryQueryStats stats;
    const int count =
        std::min<int>(probe_count, static_cast<int>(db.movable_ids.size()));
    const ea::Real delta =
        radius_bin_scale * std::min(density.bin_width(), density.bin_height());
    std::vector<int> ids = db.movable_ids;
    std::sort(ids.begin(), ids.end(),
              [&](int a, int b) {
                  const auto ra = context.wire_x[a] * context.wire_x[a] +
                                  context.wire_y[a] * context.wire_y[a] +
                                  context.exact_density_x[a] * context.exact_density_x[a] +
                                  context.exact_density_y[a] * context.exact_density_y[a];
                  const auto rb = context.wire_x[b] * context.wire_x[b] +
                                  context.wire_y[b] * context.wire_y[b] +
                                  context.exact_density_x[b] * context.exact_density_x[b] +
                                  context.exact_density_y[b] * context.exact_density_y[b];
                  return ra < rb;
              });
    for (int k = 0; k < count; ++k) {
        const int id = ids[k];
        const ea::Node& node = db.nodes[id];
        for (int axis = 0; axis < 2; ++axis) {
            ea::Real best_slope = 0.0;
            for (const int sign : {-1, 1}) {
                const ea::Real original = axis == 0 ? node.x : node.y;
                const ea::Real trial_coord = axis == 0
                    ? std::clamp(node.x + sign * delta,
                                 db.xl + 0.5 * node.width,
                                 db.xh - 0.5 * node.width)
                    : std::clamp(node.y + sign * delta,
                                 db.yl + 0.5 * node.height,
                                 db.yh - 0.5 * node.height);
                const ea::Real actual_delta = std::abs(trial_coord - original);
                if (actual_delta <= 1e-12) continue;
                const ea::Real new_x = axis == 0 ? trial_coord : node.x;
                const ea::Real new_y = axis == 0 ? node.y : trial_coord;
                const ea::DensityMove move =
                    density.evaluate_move(id, new_x, new_y);
                ++stats.local_density_queries;
                if (move.energy_delta < -descent_tolerance) {
                    best_slope = sign * move.energy_delta / actual_delta;
                    break;
                }
            }
            if (axis == 0) aux_x[id] = best_slope;
            else aux_y[id] = best_slope;
        }
    }
    return stats;
}

}  // namespace nsgp::modules
