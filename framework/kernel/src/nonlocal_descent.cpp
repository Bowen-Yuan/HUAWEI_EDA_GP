#include "epsilon_active/nonlocal_descent.hpp"

#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ea {
namespace {
Real rms(const std::vector<Real>& x, const std::vector<Real>& y) {
    Real sum = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) sum += x[i] * x[i] + y[i] * y[i];
    return std::sqrt(sum / std::max<std::size_t>(1, 2 * x.size()));
}
void normalize(std::vector<Real>& x, std::vector<Real>& y) {
    const Real scale = rms(x, y);
    if (!std::isfinite(scale)) throw std::runtime_error("non-finite direction");
    if (scale <= 1.0e-18) return;
    for (std::size_t i = 0; i < x.size(); ++i) { x[i] /= scale; y[i] /= scale; }
}
Real active_rms(const std::vector<Real>& x, const std::vector<Real>& y,
                Real* fraction) {
    Real sum = 0.0; std::size_t active = 0;
    for (std::size_t i=0;i<x.size();++i) { const Real q=x[i]*x[i]+y[i]*y[i];
        if (q>1.0e-24) { sum+=q; ++active; } }
    if (fraction) *fraction=static_cast<Real>(active)/std::max<std::size_t>(1,x.size());
    return active ? std::sqrt(sum/(2.0*active)) : 0.0;
}
void clamp(Database& db) {
    for (int id : db.movable_ids) {
        Node& n = db.nodes[id];
        n.x = std::clamp(n.x, db.xl + .5 * n.width, db.xh - .5 * n.width);
        n.y = std::clamp(n.y, db.yl + .5 * n.height, db.yh - .5 * n.height);
    }
}
}

NonlocalDescentStats run_nonlocal_descent(Database& db, int bins_x, int bins_y,
    Real target_density, const NonlocalDescentConfig& c,
    const AuxiliaryGradient& auxiliary, const NonlocalTelemetry& telemetry) {
    if (c.iterations < 1 || c.learning_rate <= 0.0 || c.lambda_minimum <= 0.0 ||
        c.lambda_minimum > c.lambda_maximum || c.lambda_update_interval < 1)
        throw std::invalid_argument("invalid nonlocal pure-descent configuration");
    ExactHpwl hpwl(db);
    ExactOverlapDensity density(db, bins_x, bins_y, target_density);
    auto optimizer = make_optimizer(c.optimizer, c.beta1, c.beta2, .7, c.numerical_epsilon);
    const std::size_t n = db.movable_ids.size();
    optimizer->reset(2 * n);
    std::vector<Real> wx, wy, dx, dy, ax, ay, gradient(2 * n), delta;
    NonlocalDescentStats result;
    LambdaConfig lambda_config;
    lambda_config.policy=LambdaPolicy::Trajectory;
    lambda_config.density_weight_scale=1.0;
    lambda_config.initial_effective=c.lambda_initial;
    lambda_config.stop_overflow=c.stop_overflow_percent*0.01;
    lambda_config.trajectory_start_overflow=c.start_overflow_percent*0.01;
    lambda_config.update_interval=c.lambda_update_interval;
    lambda_config.trajectory_horizon=c.iterations;
    lambda_config.control_min=c.lambda_minimum;
    lambda_config.control_max=c.lambda_maximum;
    LambdaController controller(lambda_config); bool initialized=false;
    const Real die_scale = std::min(db.xh - db.xl, db.yh - db.yl);
    for (int iteration = 0; iteration < c.iterations; ++iteration) {
        const DensityMetrics metrics = density.evaluate(0.0, 1.0, &dx, &dy);
        hpwl.evaluate(c.hpwl_epsilon * std::min(density.bin_width(), density.bin_height()),
                      4.0, 100, &wx, &wy);
        hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        if (!initialized) { controller.initialize(db,wx,wy,dx,dy,0.0,metrics.overflow); initialized=true; }
        else controller.update(iteration,0.0,metrics.overflow);
        const Real lambda=controller.effective();
        ax.assign(db.nodes.size(), 0.0); ay.assign(db.nodes.size(), 0.0);
        AuxiliaryContext aux_context{iteration,lambda,wx,wy,dx,dy,0.0,metrics};
        auxiliary(db, density, aux_context, ax, ay);
        Real aux_fraction=0.0; const Real aux_rms = active_rms(ax, ay, &aux_fraction);
        if (aux_rms > 1.0e-18) ++result.nonzero_auxiliary_steps;
        normalize(wx, wy); normalize(dx, dy);
        if (aux_rms>1.0e-18) for(std::size_t i=0;i<ax.size();++i){ax[i]/=aux_rms;ay[i]/=aux_rms;}
        for (std::size_t m = 0; m < n; ++m) {
            const int id = db.movable_ids[m];
            gradient[m] = wx[id] + lambda * (c.exact_density_weight * dx[id] +
                c.auxiliary_density_weight * ax[id]);
            gradient[n + m] = wy[id] + lambda * (c.exact_density_weight * dy[id] +
                c.auxiliary_density_weight * ay[id]);
            if (!std::isfinite(gradient[m]) || !std::isfinite(gradient[n + m]))
                throw std::runtime_error("non-finite combined direction");
        }
        optimizer->compute_delta(gradient, c.learning_rate * die_scale,
                                 c.maximum_delta * c.learning_rate * die_scale, delta);
        for (std::size_t m = 0; m < n; ++m) {
            Node& node = db.nodes[db.movable_ids[m]];
            node.x -= delta[m]; node.y -= delta[n + m];
        }
        clamp(db);
        if (telemetry) {
            const DensityMetrics after = density.evaluate(0.0, 1.0, nullptr, nullptr);
            const Real after_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
            telemetry(iteration, after_hpwl, after, lambda, rms(wx, wy), rms(dx, dy), aux_rms);
            result.objective_evaluations += 2;
        }
        result.maximum_overflow = std::max(result.maximum_overflow, metrics.overflow);
        result.objective_evaluations += 2;
    }
    return result;
}
}  // namespace ea
