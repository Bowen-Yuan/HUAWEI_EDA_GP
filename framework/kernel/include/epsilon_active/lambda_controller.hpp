#pragma once

#include "epsilon_active/types.hpp"

#include <string>
#include <vector>

namespace ea {

enum class LambdaPolicy { Dreamplace, Trajectory, Ratio };

LambdaPolicy parse_lambda_policy(const std::string& name);
const char* lambda_policy_name(LambdaPolicy policy) noexcept;

struct LambdaConfig {
    LambdaPolicy policy = LambdaPolicy::Dreamplace;
    Real density_weight_scale = 8.0e-5;
    Real reference_hpwl_delta = 3.5e5;
    Real stop_overflow = 0.07;
    Real hpwl_baseline = 0.0;
    Real overflow_baseline = 0.07;
    // A non-negative value is the soft trajectory reference start.  The
    // historical default (-1) retains initial-overflow behavior.
    Real trajectory_start_overflow = -1.0;
    int update_interval = 25;
    int trajectory_horizon = 550;
    Real control_min = 1.0e-6;
    Real control_max = 1.0e12;
    // Optional staged-continuation value. Zero keeps the normal fresh
    // initialization; a positive value restores the previous effective
    // multiplier after the new level recomputes its gradient scale.
    Real initial_effective = 0.0;
    Real initial_control = 1.0;
};

class LambdaController {
public:
    explicit LambdaController(LambdaConfig config);

    void initialize(const Database& db, const std::vector<Real>& wire_x,
                    const std::vector<Real>& wire_y,
                    const std::vector<Real>& density_x,
                    const std::vector<Real>& density_y, Real hpwl,
                    Real overflow);
    void update(int density_step, Real hpwl, Real overflow);

    Real base() const noexcept { return base_; }
    Real control() const noexcept { return control_; }
    Real effective() const noexcept { return base_ * control_; }

private:
    LambdaConfig config_;
    Real base_ = 0.0;
    Real control_ = 1.0;
    Real previous_hpwl_ = 0.0;
    Real previous_overflow_ = 0.0;
    Real initial_overflow_ = 0.0;
    Real integral_error_ = 0.0;
};

}  // namespace ea
