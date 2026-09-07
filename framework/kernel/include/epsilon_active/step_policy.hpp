#pragma once

#include "epsilon_active/types.hpp"

#include <memory>
#include <string>

namespace ea {

// Lightweight step controllers sit between the optimizer and the search
// loop.  They only decide the learning rate and the per-step maximum delta;
// acceptance and backtracking stay in the caller so every optimizer is
// compared under the same exact contract.
enum class StepPolicyKind {
    Constant,
    CosineDecay,
    TrustRadius
};

StepPolicyKind parse_step_policy(const std::string& name);
const char* step_policy_name(StepPolicyKind kind) noexcept;

struct StepDecision {
    Real learning_rate = 0.0;
    Real maximum_delta = 0.0;
};

struct StepObservation {
    bool accepted = false;
    int backtracks = 0;
};

// One configuration struct keeps make_step_controller trivial.  Unused
// fields are ignored by each policy; validate() checks every field that any
// policy reads.
struct StepControllerConfig {
    Real learning_rate = 0.02;
    Real maximum_delta_bins = 0.5;
    Real learning_rate_min_ratio = 0.10;
    Real trust_radius_bins = 1.0;
    Real radius_min_bins = 0.02;
    Real radius_max_bins = 2.0;
    Real grow_factor = 1.25;
    Real shrink_factor = 0.5;
    int grow_after_accepts = 3;

    void validate() const;
};

class StepController {
public:
    virtual ~StepController() = default;

    virtual void reset(int total_iterations, Real bin_size) = 0;
    virtual StepDecision propose(int iteration) = 0;
    virtual void observe(const StepObservation& observation) = 0;

    // Clears accept/reject streaks without changing the radius or decay
    // schedule.  Used when the optimizer state is reset.
    virtual void reset_streaks() {}
};

std::unique_ptr<StepController> make_step_controller(
    StepPolicyKind kind, const StepControllerConfig& config);

}  // namespace ea
