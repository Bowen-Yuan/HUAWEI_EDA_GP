#pragma once

#include "microkernel.hpp"
#include "epsilon_active/nonlocal_descent.hpp"

namespace nsgp::nonlocal {
inline ea::NonlocalDescentConfig config_from_json(const Json& j) {
    ea::NonlocalDescentConfig c;
    c.iterations = j.value("iterations", c.iterations);
    if (j.contains("optimizer")) { const auto& o=j.at("optimizer");
        c.optimizer=ea::parse_optimizer(o.value("name",std::string("amsgrad")));
        c.learning_rate=o.value("learning_rate",c.learning_rate);
        c.beta1=o.value("beta1",c.beta1); c.beta2=o.value("beta2",c.beta2);
        c.numerical_epsilon=o.value("numerical_epsilon",c.numerical_epsilon); }
    if (j.contains("direction")) { const auto& d=j.at("direction");
        c.hpwl_epsilon=d.value("wire_epsilon_bin_scale",c.hpwl_epsilon);
        c.exact_density_weight=d.value("exact_density_weight",c.exact_density_weight);
        c.auxiliary_density_weight=d.value("aux_density_weight",c.auxiliary_density_weight); }
    if (j.contains("lambda")) { const auto& l=j.at("lambda");
        const auto policy=l.value("policy",std::string("trajectory"));
        if (policy != "trajectory") throw std::invalid_argument("nonlocal modules require lambda.policy=trajectory");
        c.lambda_initial=l.value("initial",c.lambda_initial); c.lambda_minimum=l.value("minimum",c.lambda_minimum);
        c.lambda_maximum=l.value("maximum",c.lambda_maximum);
        c.start_overflow_percent=l.value("start_overflow_percent",c.start_overflow_percent);
        c.stop_overflow_percent=l.value("stop_overflow_percent",c.stop_overflow_percent);
        c.lambda_update_interval=l.value("update_interval",c.lambda_update_interval);
        c.lambda_kp=l.value("kp",c.lambda_kp); c.lambda_ki=l.value("ki",c.lambda_ki); c.lambda_kd=l.value("kd",c.lambda_kd); }
    return c;
}
inline StageStats run(StageContext& x, const Json& j, const ea::AuxiliaryDirection& aux) {
    const auto c=config_from_json(j);
    const auto s=ea::run_nonlocal_descent(x.db,x.density.bins_x,x.density.bins_y,x.density.target_density,c,aux,
        [&](int iteration, ea::Real hpwl, const ea::DensityMetrics& d, ea::Real lambda,
            ea::Real wire_rms, ea::Real exact_rms, ea::Real aux_rms) {
            if (x.log_iteration) x.log_iteration(iteration, "nonlocal", Json{
                {"hpwl",hpwl},{"overflow_percent",d.overflow*100.0},{"density_energy",d.energy},
                {"max_density",d.max_density},{"lambda",lambda},{"wire_grad_rms",wire_rms},
                {"exact_density_grad_rms",exact_rms},{"aux_density_grad_rms",aux_rms},
                {"objective_evaluations",4}}); });
    return {c.iterations,c.iterations,0,s.objective_evaluations};
}
}  // namespace nsgp::nonlocal
