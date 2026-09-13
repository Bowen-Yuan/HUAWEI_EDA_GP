#pragma once

#include "epsilon_active/nonlocal_descent.hpp"

namespace nsgp::modules {

// Exact finite-radius secant probe.  Only strictly descending directions
// (energy_delta < -descent_tolerance) contribute; the secant denominator is
// the actual (post-clamp) displacement.  Probes never commit: the shared
// occupancy must correspond to the live layout before the call.
ea::AuxiliaryQueryStats finite_radius_probe(
    const ea::Database& db, ea::ExactOverlapDensity& density,
    const ea::AuxiliaryContext& context,
    int probe_count, ea::Real radius_bin_scale, ea::Real descent_tolerance,
    std::vector<ea::Real>& aux_x, std::vector<ea::Real>& aux_y);

}  // namespace nsgp::modules
