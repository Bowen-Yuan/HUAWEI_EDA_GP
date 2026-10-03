#include "epsilon_active/batch_acceptance.hpp"
#include "epsilon_active/bisection.hpp"
#include "epsilon_active/coarse_flow.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/density_coordinate.hpp"
#include "epsilon_active/hpwl.hpp"
#include "epsilon_active/optimizer.hpp"
#include "epsilon_active/recovery.hpp"
#include "epsilon_active/swap_recovery.hpp"
#include "epsilon_active/transport.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(std::string(message) + ": actual=" +
                                 std::to_string(actual) +
                                 " expected=" + std::to_string(expected));
    }
}

ea::Database make_database() {
    ea::Database db;
    db.xl = 0.0;
    db.yl = 0.0;
    db.xh = 20.0;
    db.yh = 10.0;
    db.nodes = {
        {0, "a", 5.0, 5.0, 10.0, 10.0, false, false, "N"},
        {1, "b", 15.0, 5.0, 2.0, 2.0, false, false, "N"},
        {2, "macro", 5.0, 5.0, 10.0, 10.0, true, false, "N"},
        {3, "port", 5.0, 5.0, 10.0, 10.0, true, true, "N"},
    };
    db.movable_ids = {0, 1};
    db.fixed_ids = {2, 3};
    db.movable_area = 104.0;
    db.pins = {{0, 0.0, 0.0}, {1, 0.0, 0.0}};
    db.nets = {{0, "n0", 0, 2, 1.0}};
    db.node_pin_count = {1, 1, 0, 0};
    db.node_pin_offsets = {0, 1, 2, 2, 2};
    db.node_pin_indices = {0, 1};
    return db;
}

void test_hpwl() {
    ea::Database db = make_database();
    ea::ExactHpwl oracle(db);
    std::vector<ea::Real> gx, gy;
    const double value = oracle.evaluate(0.0, 4.0, 100, &gx, &gy);
    near(value, 10.0, 1.0e-12, "exact HPWL");
    near(gx[0], -1.0, 1.0e-12, "left pin subgradient");
    near(gx[1], 1.0, 1.0e-12, "right pin subgradient");
    near(gy[0] + gy[1], 0.0, 1.0e-12, "tied y subgradient balance");

    const double active_value = oracle.evaluate(20.0, 4.0, 100, &gx, &gy);
    near(active_value, value, 1.0e-12, "epsilon must not change HPWL value");
    near(gx[0] + gx[1], 0.0, 1.0e-12, "active HPWL balance");
}

void test_fixed_macro_density() {
    ea::Database db = make_database();
    ea::ExactOverlapDensity density(db, 2, 1, 1.0);
    std::vector<ea::Real> gx, gy;
    const ea::DensityMetrics blocked = density.evaluate(0.0, 4.0, &gx, &gy);
    near(blocked.overflow, 100.0 / 104.0, 1.0e-12,
         "fixed macro must consume bin capacity");
    near(density.fixed_occupancy()[0], 100.0, 1.0e-12,
         "fixed macro occupancy");
    near(density.fixed_occupancy()[1], 0.0, 1.0e-12,
         "FIXED_NI must not consume capacity");

    db.nodes[0].x = 15.0;
    const ea::DensityMetrics clear = density.evaluate(0.0, 4.0, &gx, &gy);
    near(clear.overflow, 4.0 / 104.0, 1.0e-12,
         "movable clear of fixed macro");
}

void test_overlap_direction_and_exact_value() {
    ea::Database db = make_database();
    db.nodes[0].x = 9.0;
    db.nodes[0].width = 4.0;
    db.movable_area = 44.0;
    ea::ExactOverlapDensity density(db, 2, 1, 1.0);
    std::vector<ea::Real> exact_x, exact_y, active_x, active_y;
    const auto exact = density.evaluate(0.0, 4.0, &exact_x, &exact_y);
    const auto active = density.evaluate(3.0, 4.0, &active_x, &active_y);
    const auto linear_active = density.evaluate(3.0, 1.0, &active_x, &active_y);
    near(active.energy, exact.energy, 1.0e-12,
         "epsilon must not change overlap energy");
    near(active.overflow, exact.overflow, 1.0e-12,
         "epsilon must not change overflow");
    near(linear_active.energy, exact.energy, 1.0e-12,
         "density active power must not change overlap energy");
    near(linear_active.overflow, exact.overflow, 1.0e-12,
         "density active power must not change overflow");
    require(exact_x[0] < 0.0, "overfull fixed-macro bin must push cell right");

    std::vector<ea::Real> weighted_x, weighted_y;
    const std::vector<ea::Real> uniform_prices(2, 1.0);
    const auto weighted = density.evaluate_with_prices(
        uniform_prices, 0.0, 4.0, &weighted_x, &weighted_y);
    near(weighted.overflow, exact.overflow, 1.0e-12,
         "uniform regional prices must preserve exact overflow");
    near(weighted_x[0], exact_x[0], 1.0e-12,
         "uniform regional prices must preserve x direction");
    near(weighted_y[0], exact_y[0], 1.0e-12,
         "uniform regional prices must preserve y direction");
}

void test_bin_field_prolongation() {
    const std::vector<double> coarse{1.0, 2.0, 3.0, 4.0};
    const std::vector<double> fine =
        ea::prolongate_bin_field(coarse, 2, 2, 4, 4);
    require(fine.size() == 16, "prolongated field must have target size");
    near(fine[0], 1.0, 1.0e-12,
         "top-left child must inherit its exact parent price");
    near(fine[3], 2.0, 1.0e-12,
         "top-right child must inherit its exact parent price");
    near(fine[12], 3.0, 1.0e-12,
         "bottom-left child must inherit its exact parent price");
    near(fine[15], 4.0, 1.0e-12,
         "bottom-right child must inherit its exact parent price");
}

void test_optimizer_factory() {
    for (const char* name : {"adam", "amsgrad", "adagrad", "heavy-ball", "sgd"}) {
        auto optimizer = ea::make_optimizer(
            ea::parse_optimizer(name), 0.9, 0.99, 0.7, 1.0e-8);
        std::vector<ea::Real> delta;
        optimizer->compute_delta({1.0, -1.0}, 0.1, 0.2, delta);
        require(delta.size() == 2, "optimizer output dimension");
        require(delta[0] > 0.0 && delta[1] < 0.0,
                "optimizer descent direction signs");
    }
}

void test_exact_batch_backtracking() {
    std::vector<ea::Real> state = {0.0};
    ea::BatchAcceptanceConfig config;
    config.enabled = true;
    config.max_trials = 4;
    config.shrink = 0.5;
    config.overflow_cap = 0.07;
    const auto apply = [&](const std::vector<ea::Real>& candidate) {
        state = candidate;
    };
    const auto evaluate = [&]() {
        const double hpwl = (state[0] - 0.4) * (state[0] - 0.4);
        const double overflow = state[0] > 0.75 ? 0.10 : 0.05;
        return std::make_pair(hpwl, overflow);
    };
    const auto accepted = ea::accept_exact_batch(
        state, {-1.0}, 0.16, 0.05, config, apply, evaluate);
    require(accepted.accepted, "smaller exact batch must be accepted");
    near(accepted.scale, 0.5, 1.0e-12, "accepted half step");
    near(state[0], 0.5, 1.0e-12, "accepted position");

    const std::vector<ea::Real> saved = state;
    const auto rejected = ea::accept_exact_batch(
        saved, {-1.0}, accepted.hpwl, accepted.overflow, config, apply,
        [&]() { return std::make_pair(1.0, 0.10); });
    require(!rejected.accepted, "infeasible batches must be rejected");
    require(state == saved, "rejected batch must restore positions exactly");

    ea::BatchAcceptanceConfig infeasible = config;
    infeasible.allow_infeasible = true;
    infeasible.max_infeasible_hpwl_increase = 0.10;
    state = {0.0};
    const auto exploratory = ea::accept_exact_batch(
        state, {-1.0}, 0.16, 0.20, infeasible, apply,
        [&]() { return std::make_pair(0.17, state[0] > 0.0 ? 0.10 : 0.20); });
    require(exploratory.accepted,
            "infeasible exact batch may accept strict overflow descent");
    near(state[0], 1.0, 1.0e-12,
         "infeasible exact batch must commit the audited candidate");
}

void test_exact_capacity_transport() {
    ea::Database db = make_database();
    ea::ExactOverlapDensity density(db, 2, 1, 1.0);
    ea::TransportConfig config;
    config.rounds = 2;
    config.max_moves = 4;
    config.max_source_bins = 2;
    config.candidate_lookahead = 2;
    const double fixed_x = db.nodes[2].x;
    const ea::TransportStats stats = ea::transport_excess_to_capacity(db, density, config);
    require(stats.moves > 0, "transport must move an exact-overflow contributor");
    require(stats.final_overflow < stats.initial_overflow,
            "every transport result must reduce exact overflow");
    near(db.nodes[2].x, fixed_x, 1.0e-12, "transport must not move fixed macros");

    ea::Database weighted_db = make_database();
    ea::ExactOverlapDensity weighted_density(weighted_db, 2, 1, 1.0);
    config.group_size = 2;
    config.weighted_grouping = true;
    const ea::TransportStats weighted_stats =
        ea::transport_excess_to_capacity(weighted_db, weighted_density, config);
    require(weighted_stats.moves > 0 &&
            weighted_stats.final_overflow < weighted_stats.initial_overflow,
            "weighted groups must preserve exact-overflow transport descent");

    ea::Database offset_db = make_database();
    offset_db.nodes[0].x = 3.0;
    offset_db.nodes[0].y = 3.0;
    offset_db.nodes[0].width = 2.0;
    offset_db.nodes[0].height = 2.0;
    offset_db.movable_area = 8.0;
    ea::ExactOverlapDensity offset_density(offset_db, 2, 1, 1.0);
    config.group_size = 1;
    config.weighted_grouping = false;
    config.preserve_bin_offset = true;
    const ea::TransportStats offset_stats =
        ea::transport_excess_to_capacity(offset_db, offset_density, config);
    require(offset_stats.moves > 0 &&
            offset_stats.final_overflow < offset_stats.initial_overflow,
            "offset-preserving transport must retain exact-overflow descent");
    near(offset_db.nodes[0].x, 13.0, 1.0e-12,
         "transport must preserve x offset from bin center");
    near(offset_db.nodes[0].y, 3.0, 1.0e-12,
         "transport must preserve y offset from bin center");

    ea::Database rigid_db = make_database();
    rigid_db.nodes[0].x = 3.0;
    rigid_db.nodes[0].y = 3.0;
    rigid_db.nodes[0].width = 2.0;
    rigid_db.nodes[0].height = 2.0;
    rigid_db.nodes[1].x = 7.0;
    rigid_db.nodes[1].y = 7.0;
    rigid_db.movable_area = 8.0;
    ea::ExactOverlapDensity rigid_density(rigid_db, 2, 1, 1.0);
    config.group_size = 2;
    config.preserve_bin_offset = true;
    config.rigid_group_displacement = true;
    const double relative_x = rigid_db.nodes[1].x - rigid_db.nodes[0].x;
    const double relative_y = rigid_db.nodes[1].y - rigid_db.nodes[0].y;
    const ea::TransportStats rigid_stats =
        ea::transport_excess_to_capacity(rigid_db, rigid_density, config);
    require(rigid_stats.moves >= 2 &&
            rigid_stats.final_overflow < rigid_stats.initial_overflow,
            "rigid group transport must move both nodes with exact descent");
    near(rigid_db.nodes[1].x - rigid_db.nodes[0].x, relative_x, 1.0e-12,
         "rigid transport must preserve relative x");
    near(rigid_db.nodes[1].y - rigid_db.nodes[0].y, relative_y, 1.0e-12,
         "rigid transport must preserve relative y");

    ea::Database atomic_db;
    atomic_db.xl = 0.0;
    atomic_db.yl = 0.0;
    atomic_db.xh = 30.0;
    atomic_db.yh = 10.0;
    atomic_db.nodes = {
        {0, "a", 3.0, 5.0, 2.0, 2.0, false, false, "N"},
        {1, "b", 13.0, 5.0, 2.0, 2.0, false, false, "N"},
        {2, "left_macro", 5.0, 5.0, 10.0, 10.0, true, false, "N"},
        {3, "right_macro", 25.0, 5.0, 9.8, 10.0, true, false, "N"},
    };
    atomic_db.movable_ids = {0, 1};
    atomic_db.fixed_ids = {2, 3};
    atomic_db.movable_area = 8.0;
    atomic_db.pins = {{0, 0.0, 0.0}, {1, 0.0, 0.0}};
    atomic_db.nets = {{0, "pair", 0, 2, 1.0}};
    atomic_db.node_pin_count = {1, 1, 0, 0};
    atomic_db.node_pin_offsets = {0, 1, 2, 2, 2};
    atomic_db.node_pin_indices = {0, 1};
    ea::ExactOverlapDensity atomic_density(atomic_db, 3, 1, 1.0);
    atomic_density.evaluate(0.0, 1.0, nullptr, nullptr);
    const ea::DensityMove blocked_member = atomic_density.evaluate_move(1, 23.0, 5.0);
    require(blocked_member.overflow_area_delta > 0.0,
            "second atomic member must be blocked by member-wise descent");
    config.rounds = 1;
    config.max_source_bins = 3;
    config.group_size = 2;
    config.rigid_group_displacement = false;
    config.atomic_groups = true;
    const double left_macro_x = atomic_db.nodes[2].x;
    const double right_macro_x = atomic_db.nodes[3].x;
    const ea::TransportStats atomic_stats =
        ea::transport_excess_to_capacity(atomic_db, atomic_density, config);
    require(atomic_stats.atomic_attempts == 1 && atomic_stats.atomic_groups == 1,
            "atomic group transaction must be attempted and accepted once");
    require(atomic_stats.moves == 2 &&
            atomic_stats.final_overflow < atomic_stats.initial_overflow,
            "atomic group must move both members under aggregate exact descent");
    near(atomic_db.nodes[0].x, 13.0, 1.0e-12,
         "first atomic member displacement");
    near(atomic_db.nodes[1].x, 23.0, 1.0e-12,
         "positive-delta member must move with its atomic group");
    near(atomic_db.nodes[2].x, left_macro_x, 1.0e-12,
         "atomic transport must not move the left fixed macro");
    near(atomic_db.nodes[3].x, right_macro_x, 1.0e-12,
         "atomic transport must not move the right fixed macro");

    ea::Database rejected_db = atomic_db;
    rejected_db.nodes[0].x = 3.0;
    rejected_db.nodes[1].x = 13.0;
    rejected_db.nodes[3].width = 10.0;
    ea::ExactOverlapDensity rejected_density(rejected_db, 3, 1, 1.0);
    const ea::TransportStats rejected_stats =
        ea::transport_excess_to_capacity(rejected_db, rejected_density, config);
    require(rejected_stats.atomic_attempts == 1 && rejected_stats.atomic_groups == 0,
            "zero-delta atomic group must be rejected");
    require(rejected_stats.moves == 0,
            "rejected atomic group must not count member moves");
    near(rejected_db.nodes[0].x, 3.0, 1.0e-12,
         "rejected atomic group must restore first member");
    near(rejected_db.nodes[1].x, 13.0, 1.0e-12,
         "rejected atomic group must restore positive-delta member");
    near(rejected_stats.final_overflow, rejected_stats.initial_overflow, 1.0e-12,
         "rejected atomic group must restore exact overflow");

    ea::Database fallback_db = atomic_db;
    fallback_db.nodes[0].x = 3.0;
    fallback_db.nodes[1].x = 13.0;
    ea::ExactOverlapDensity fallback_density(fallback_db, 3, 1, 1.0);
    config.atomic_min_gain_ratio = 0.75;
    config.atomic_fallback_individual = true;
    const ea::TransportStats fallback_stats =
        ea::transport_excess_to_capacity(fallback_db, fallback_density, config);
    require(fallback_stats.atomic_attempts == 1 &&
            fallback_stats.atomic_groups == 0 &&
            fallback_stats.atomic_fallbacks == 1,
            "density-inefficient atomic group must fall back once");
    require(fallback_stats.moves == 1 &&
            fallback_stats.final_overflow < fallback_stats.initial_overflow,
            "individual fallback must retain the first member's exact descent");
    near(fallback_db.nodes[0].x, 13.0, 1.0e-12,
         "individual fallback must apply the first member move");
    near(fallback_db.nodes[1].x, 13.0, 1.0e-12,
         "individual fallback must restore the other atomic member");

    ea::Database active_db;
    active_db.xl = 0.0;
    active_db.yl = 0.0;
    active_db.xh = 30.0;
    active_db.yh = 10.0;
    active_db.nodes = {
        {0, "a", 3.0, 5.0, 2.0, 2.0, false, false, "N"},
        {1, "c", 7.0, 5.0, 2.0, 2.0, false, false, "N"},
        {2, "remote", 13.0, 5.0, 2.0, 2.0, false, false, "N"},
        {3, "left_macro", 5.0, 5.0, 10.0, 10.0, true, false, "N"},
    };
    active_db.movable_ids = {0, 1, 2};
    active_db.fixed_ids = {3};
    active_db.movable_area = 12.0;
    active_db.pins = {{0, 0.0, 0.0}, {1, 0.0, 0.0}, {2, 0.0, 0.0}};
    active_db.nets = {{0, "triple", 0, 3, 1.0}};
    active_db.node_pin_count = {1, 1, 1, 0};
    active_db.node_pin_offsets = {0, 1, 2, 3, 3};
    active_db.node_pin_indices = {0, 1, 2};
    ea::ExactOverlapDensity active_density(active_db, 3, 1, 1.0);
    config.group_size = 3;
    config.atomic_min_gain_ratio = 0.0;
    config.atomic_fallback_individual = false;
    config.atomic_source_active = true;
    const ea::TransportStats active_stats =
        ea::transport_excess_to_capacity(active_db, active_density, config);
    require(active_stats.atomic_attempts == 1 && active_stats.atomic_groups == 1,
            "source-active subgroup must transact once");
    require(active_stats.moves == 2 &&
            active_stats.final_overflow < active_stats.initial_overflow,
            "only two exact source contributors must move atomically");
    near(active_db.nodes[0].x, 13.0, 1.0e-12,
         "first source-active member displacement");
    near(active_db.nodes[1].x, 17.0, 1.0e-12,
         "second source-active member displacement");
    near(active_db.nodes[2].x, 13.0, 1.0e-12,
         "remote group member must remain eligible and unmoved");
    near(active_db.nodes[3].x, 5.0, 1.0e-12,
         "source-active atomic transport must not move fixed macros");

    ea::Database scheduled_db;
    scheduled_db.xl = 0.0;
    scheduled_db.yl = 0.0;
    scheduled_db.xh = 40.0;
    scheduled_db.yh = 10.0;
    scheduled_db.nodes = {
        {0, "a", 3.0, 5.0, 2.0, 2.0, false, false, "N"},
        {1, "b", 7.0, 5.0, 2.0, 2.0, false, false, "N"},
        {2, "source_macro", 5.0, 5.0, 10.0, 10.0, true, false, "N"},
        {3, "destination_macro", 15.0, 5.0, 9.6, 10.0, true, false, "N"},
    };
    scheduled_db.movable_ids = {0, 1};
    scheduled_db.fixed_ids = {2, 3};
    scheduled_db.movable_area = 8.0;
    scheduled_db.pins = {{0, 0.0, 0.0}, {1, 0.0, 0.0}};
    scheduled_db.nets = {{0, "pair", 0, 2, 1.0}};
    scheduled_db.node_pin_count = {1, 1, 0, 0};
    scheduled_db.node_pin_offsets = {0, 1, 2, 2, 2};
    scheduled_db.node_pin_indices = {0, 1};
    ea::ExactOverlapDensity scheduled_density(scheduled_db, 4, 1, 1.0);
    config.rounds = 2;
    config.max_moves = 4;
    config.max_source_bins = 4;
    config.group_size = 2;
    config.atomic_source_active = true;
    config.atomic_rounds = 1;
    const ea::TransportStats scheduled_stats =
        ea::transport_excess_to_capacity(scheduled_db, scheduled_density, config);
    require(scheduled_stats.atomic_attempts == 1 && scheduled_stats.atomic_groups == 1,
            "only the atomic prefix round may transact the source-active group");
    require(scheduled_stats.moves == 3,
            "the post-prefix round must finish with one independent exact move");
    near(scheduled_stats.final_overflow, 0.0, 1.0e-12,
         "the scheduled transport must remove exact overflow");
    near(scheduled_db.nodes[0].x, 23.0, 1.0e-12,
         "post-prefix independent move must use offset-preserving transport");
    near(scheduled_db.nodes[1].x, 17.0, 1.0e-12,
         "post-prefix cleanup must leave the second group member in place");

    ea::Database guarded_db = scheduled_db;
    guarded_db.nodes[0].x = 3.0;
    guarded_db.nodes[1].x = 7.0;
    ea::ExactOverlapDensity guarded_density(guarded_db, 4, 1, 1.0);
    config.rounds = 1;
    config.max_moves = 1;
    config.atomic_rounds = -1;
    config.atomic_componentwise_capacity = true;
    config.atomic_fallback_individual = true;
    const ea::TransportStats guarded_stats =
        ea::transport_excess_to_capacity(guarded_db, guarded_density, config);
    require(guarded_stats.atomic_attempts == 1 && guarded_stats.atomic_groups == 0 &&
            guarded_stats.atomic_capacity_rejections == 1,
            "atomic destination spill must fail the componentwise capacity guard");
    require(guarded_stats.atomic_fallbacks == 1 && guarded_stats.moves == 1,
            "capacity rejection must use the exact individual fallback");
    near(guarded_db.nodes[0].x, 13.0, 1.0e-12,
         "guarded fallback must commit the audited first-member move");
    near(guarded_db.nodes[1].x, 7.0, 1.0e-12,
         "guarded fallback must restore the remaining atomic member");

    ea::Database bidding_db = scheduled_db;
    bidding_db.nodes[0].x = 3.0;
    bidding_db.nodes[1].x = 7.0;
    ea::ExactOverlapDensity bidding_density(bidding_db, 4, 1, 1.0);
    config.max_moves = 4;
    config.atomic_capacity_candidates = 2;
    const ea::TransportStats bidding_stats =
        ea::transport_excess_to_capacity(bidding_db, bidding_density, config);
    require(bidding_stats.atomic_attempts == 1 && bidding_stats.atomic_groups == 1,
            "multi-candidate bidding must accept a capacity-safe atomic destination");
    require(bidding_stats.atomic_bid_trials == 2 &&
            bidding_stats.atomic_feasible_bids == 1 &&
            bidding_stats.atomic_rescued_groups == 1,
            "a later exact-capacity bid must rescue the failed first candidate");
    require(bidding_stats.atomic_fallbacks == 0 && bidding_stats.moves == 2,
            "a rescued group must not use individual fallback");
    near(bidding_db.nodes[0].x, 23.0, 1.0e-12,
         "capacity bidding must commit the exact feasible group displacement");
    near(bidding_db.nodes[1].x, 27.0, 1.0e-12,
         "capacity bidding must preserve group relative geometry");

    ea::Database batch_db = scheduled_db;
    batch_db.nodes[0].x = 3.0;
    batch_db.nodes[1].x = 7.0;
    ea::ExactOverlapDensity batch_density(batch_db, 4, 1, 1.0);
    config.atomic_source_batch = true;
    const ea::TransportStats batch_stats =
        ea::transport_excess_to_capacity(batch_db, batch_density, config);
    require(batch_stats.atomic_batch_sources == 1 &&
            batch_stats.atomic_batch_groups == 1 &&
            batch_stats.atomic_batch_bids == 1 &&
            batch_stats.atomic_batch_commits == 1,
            "source batch must generate and commit the exact feasible group bid");
    require(batch_stats.atomic_groups == 1 && batch_stats.atomic_fallbacks == 0 &&
            batch_stats.moves == 2,
            "source-batch commit must replace individual fallback");
    near(batch_db.nodes[0].x, 23.0, 1.0e-12,
         "source batch must commit the capacity-safe destination");
    near(batch_db.nodes[1].x, 27.0, 1.0e-12,
         "source batch must retain the exact rigid displacement");

    ea::Database identity_db;
    identity_db.xl = 0.0;
    identity_db.yl = 0.0;
    identity_db.xh = 30.0;
    identity_db.yh = 10.0;
    identity_db.nodes = {
        {0, "member_a", 3.0, 5.0, 2.0, 2.0, false, false, "N"},
        {1, "member_b", 7.0, 5.0, 2.0, 2.0, false, false, "N"},
        {2, "donor_a", 23.0, 5.0, 2.0, 2.0, false, false, "N"},
        {3, "donor_b", 27.0, 5.0, 2.0, 2.0, false, false, "N"},
        {4, "source_macro", 5.0, 5.0, 10.0, 10.0, true, false, "N"},
        {5, "right_anchor", 25.0, 5.0, 1.0, 1.0, true, false, "N"},
    };
    identity_db.movable_ids = {0, 1, 2, 3};
    identity_db.fixed_ids = {4, 5};
    identity_db.movable_area = 16.0;
    identity_db.pins = {
        {0, 0.0, 0.0}, {1, 0.0, 0.0},
        {0, 0.0, 0.0}, {5, 0.0, 0.0},
        {1, 0.0, 0.0}, {5, 0.0, 0.0},
        {2, 0.0, 0.0}, {4, 0.0, 0.0},
        {3, 0.0, 0.0}, {4, 0.0, 0.0},
    };
    identity_db.nets = {
        {0, "member_pair", 0, 2, 1.0},
        {1, "member_a_target", 2, 2, 5.0},
        {2, "member_b_target", 4, 2, 5.0},
        {3, "donor_a_target", 6, 2, 5.0},
        {4, "donor_b_target", 8, 2, 5.0},
    };
    identity_db.node_pin_count = {2, 2, 1, 1, 2, 2};
    identity_db.node_pin_offsets = {0, 2, 4, 5, 6, 8, 10};
    identity_db.node_pin_indices = {0, 2, 1, 4, 6, 8, 7, 9, 3, 5};
    ea::ExactOverlapDensity identity_density(identity_db, 3, 1, 1.0);
    ea::TransportConfig identity_config;
    identity_config.rounds = 0;
    identity_config.group_size = 2;
    identity_config.group_degree_limit = 32;
    identity_config.identity_exchange_passes = 1;
    identity_config.identity_exchange_candidates = 8;
    identity_config.identity_exchange_radius = 2;
    const double identity_fixed_x = identity_db.nodes[4].x;
    const ea::TransportStats identity_stats =
        ea::transport_excess_to_capacity(identity_db, identity_density,
                                         identity_config);
    require(identity_stats.identity_exchange_attempts > 0 &&
            identity_stats.identity_exchange_accepted == 1 &&
            identity_stats.identity_exchange_permuted_nodes == 4,
            "group-guided identity exchange must accept one complete permutation");
    require(identity_stats.identity_final_hpwl < identity_stats.identity_initial_hpwl,
            "identity exchange must strictly reduce exact max-minus-min HPWL");
    near(identity_stats.identity_final_overflow,
         identity_stats.identity_initial_overflow, 1.0e-12,
         "equal-shape identity exchange must preserve exact overlap overflow");
    near(identity_db.nodes[0].x, 23.0, 1.0e-12,
         "first connected identity must take the first destination anchor");
    near(identity_db.nodes[1].x, 27.0, 1.0e-12,
         "second connected identity must take the second destination anchor");
    near(identity_db.nodes[4].x, identity_fixed_x, 1.0e-12,
         "identity exchange must not move fixed macros");
}

void test_recursive_bisection_fixed_capacity() {
    ea::Database db = make_database();
    const double fixed_x = db.nodes[2].x;
    ea::ExactOverlapDensity density(db, 2, 1, 1.0);
    ea::BisectionConfig config;
    config.enabled = true;
    config.leaf_bins = 1;
    config.fm_passes = 1;
    config.leaf_hpwl_guided = true;
    config.atomic_coarsening = true;
    config.atomic_max_nodes = 2;
    config.atomic_degree_limit = 4;
    config.atomic_rounds = 1;
    config.atomic_release_depth = 1;
    const ea::BisectionStats stats =
        ea::recursive_hypergraph_bisection(db, density, config);
    require(stats.leaves > 0, "recursive bisection must create capacity leaves");
    require(stats.atomic_groups > 0,
            "atomic bisection must build deterministic supernodes");
    near(db.nodes[2].x, fixed_x, 1.0e-12,
         "recursive bisection must not move fixed macros");
    require(db.nodes[0].x >= 10.0 && db.nodes[1].x >= 10.0,
            "fixed-macro capacity must steer movable nodes to free region");
    require(std::isfinite(stats.final_hpwl) && std::isfinite(stats.final_overflow),
            "recursive bisection exact audits must be finite");
}

void test_constrained_hpwl_recovery() {
    ea::Database db = make_database();
    db.nodes[2].terminal_ni = true;
    db.nodes[0].x = 5.0;
    db.nodes[0].width = 4.0;
    db.movable_area = 28.0;
    ea::ExactOverlapDensity density(db, 2, 1, 1.0);
    ea::RecoveryConfig config;
    config.sweeps = 1;
    config.line_search_steps = 1;
    config.step_bins = 1.0;
    config.hpwl_epsilon = 0.0;
    const ea::RecoveryStats stats = ea::recover_hpwl_under_overflow(db, density, config);
    require(stats.moves > 0, "recovery must accept an exact HPWL improvement");
    require(stats.final_hpwl < stats.initial_hpwl,
            "recovery must decrease exact HPWL");
    require(stats.final_overflow <= stats.initial_overflow + 1.0e-12,
            "recovery must not increase exact overflow");
}

void test_density_neutral_swap_recovery() {
    ea::Database db = make_database();
    db.nodes[0].width = db.nodes[1].width = 2.0;
    db.nodes[0].height = db.nodes[1].height = 2.0;
    db.nodes[2].terminal_ni = true;
    db.nodes[3].x = 15.0;
    db.pins = {{0, 0.0, 0.0}, {3, 0.0, 0.0}};
    db.nets = {{0, "anchor", 0, 2, 1.0}};
    db.node_pin_count = {1, 0, 0, 1};
    db.movable_area = 8.0;
    ea::ExactOverlapDensity density(db, 2, 1, 1.0);
    ea::SwapRecoveryConfig config;
    config.sweeps = 1;
    config.radius_bins = 2;
    config.candidates = 4;
    const ea::SwapRecoveryStats stats =
        ea::recover_hpwl_with_equal_shape_swaps(db, density, config);
    require(stats.swaps > 0, "equal-shape recovery must find the anchor swap");
    require(stats.final_hpwl < stats.initial_hpwl,
            "equal-shape swap must reduce exact HPWL");
    near(stats.final_overflow, stats.initial_overflow, 1.0e-12,
         "equal-shape swap must preserve exact overflow");
}

void test_global_shape_permutation() {
    ea::Database db = make_database();
    db.yh = 20.0;
    db.nodes[0].width = db.nodes[1].width = 2.0;
    db.nodes[0].height = db.nodes[1].height = 2.0;
    db.nodes[2].terminal_ni = true;
    db.nodes[2].x = 1.0;
    db.nodes[3].x = 19.0;
    db.pins = {{0, 0.0, 0.0}, {3, 0.0, 0.0},
               {1, 0.0, 0.0}, {2, 0.0, 0.0}};
    db.nets = {{0, "right_anchor", 0, 2, 1.0},
               {1, "left_anchor", 2, 2, 1.0}};
    db.node_pin_count = {1, 1, 1, 1};
    db.movable_area = 8.0;
    ea::ExactOverlapDensity density(db, 2, 2, 1.0);
    ea::SwapRecoveryConfig config;
    config.permutation_sweeps = 1;
    const ea::SwapRecoveryStats stats =
        ea::recover_hpwl_with_equal_shape_swaps(db, density, config);
    require(stats.permutation_sweeps == 1,
            "global shape permutation must accept the crossed anchors");
    require(stats.final_hpwl < stats.initial_hpwl,
            "global shape permutation must reduce exact HPWL");
    near(stats.final_overflow, stats.initial_overflow, 1.0e-12,
         "global shape permutation must preserve exact overflow");
}

void test_net_disjoint_anchor_assignment() {
    ea::Database db;
    db.xl = 0.0;
    db.yl = 0.0;
    db.xh = 30.0;
    db.yh = 10.0;
    db.nodes = {
        {0, "a", 5.0, 5.0, 2.0, 2.0, false, false, "N"},
        {1, "b", 15.0, 5.0, 2.0, 2.0, false, false, "N"},
        {2, "c", 25.0, 5.0, 2.0, 2.0, false, false, "N"},
        {3, "ta", 15.0, 5.0, 0.0, 0.0, true, true, "N"},
        {4, "tb", 25.0, 5.0, 0.0, 0.0, true, true, "N"},
        {5, "tc", 5.0, 5.0, 0.0, 0.0, true, true, "N"},
    };
    db.movable_ids = {0, 1, 2};
    db.fixed_ids = {3, 4, 5};
    db.movable_area = 12.0;
    db.pins = {{0, 0.0, 0.0}, {3, 0.0, 0.0},
               {1, 0.0, 0.0}, {4, 0.0, 0.0},
               {2, 0.0, 0.0}, {5, 0.0, 0.0}};
    db.nets = {{0, "na", 0, 2, 1.0},
               {1, "nb", 2, 2, 1.0},
               {2, "nc", 4, 2, 1.0}};
    db.node_pin_count = {1, 1, 1, 1, 1, 1};
    db.node_pin_offsets = {0, 1, 2, 3, 4, 5, 6};
    db.node_pin_indices = {0, 2, 4, 1, 3, 5};
    ea::ExactOverlapDensity density(db, 3, 1, 1.0);
    ea::SwapRecoveryConfig config;
    config.assignment_sweeps = 1;
    config.assignment_radius_bins = 3;
    config.assignment_candidates = 3;
    config.assignment_max_bids_per_node = 128;
    const ea::SwapRecoveryStats stats =
        ea::recover_hpwl_with_equal_shape_swaps(db, density, config);
    require(stats.assignment_batches > 0,
            "net-disjoint assignment must accept an exact anchor permutation");
    require(stats.assigned_nodes == 3,
            "three-cycle anchor assignment must move every node");
    near(stats.final_hpwl, 0.0, 1.0e-12,
         "exact anchor assignment must solve the independent cycle");
    near(stats.final_overflow, stats.initial_overflow, 1.0e-12,
         "equal-shape assignment must preserve exact overflow");
}

void test_coarse_capacity_flow() {
    ea::Database db;
    db.xl = 0.0;
    db.yl = 0.0;
    db.xh = 40.0;
    db.yh = 10.0;
    db.nodes = {
        {0, "a", 5.0, 5.0, 2.0, 2.0, false, false, "N"},
        {1, "b", 15.0, 5.0, 2.0, 2.0, false, false, "N"},
        {2, "source_macro", 10.0, 5.0, 20.0, 10.0, true, false, "N"},
        {3, "target", 30.0, 5.0, 0.0, 0.0, true, true, "N"},
    };
    db.movable_ids = {0, 1};
    db.fixed_ids = {2, 3};
    db.movable_area = 8.0;
    db.pins = {{0, 0.0, 0.0}, {3, 0.0, 0.0},
               {1, 0.0, 0.0}, {3, 0.0, 0.0}};
    db.nets = {{0, "na", 0, 2, 1.0}, {1, "nb", 2, 2, 1.0}};
    db.node_pin_count = {1, 1, 0, 2};
    ea::ExactOverlapDensity density(db, 4, 1, 1.0);
    ea::CoarseFlowConfig config;
    config.enabled = true;
    config.bins_x = 2;
    config.bins_y = 1;
    config.passes = 1;
    config.exact_anchors = true;
    const double fixed_x = db.nodes[2].x;
    const ea::CoarseFlowStats stats =
        ea::coarse_capacity_flow(db, density, config);
    require(stats.flow_edges == 1 && stats.moves == 2,
            "coarse flow must realize the single source-to-sink quota");
    require(stats.exact_anchor_trials > 0 && stats.exact_anchor_rejections == 0,
            "exact-anchor flow must evaluate and accept live overlap moves");
    require(stats.final_hpwl < stats.initial_hpwl,
            "exact-HPWL realization must prefer the target-side anchors");
    require(stats.final_overflow < stats.initial_overflow,
            "coarse flow must reduce exact rectangle overlap overflow");
    near(stats.final_overflow, 0.0, 1.0e-12,
         "coarse flow must evacuate the fixed-macro source");
    near(db.nodes[2].x, fixed_x, 1.0e-12,
         "coarse flow must not move fixed macros");
}

void test_net_block_exact_recovery() {
    ea::Database db;
    db.xl = 0.0;
    db.yl = 0.0;
    db.xh = 40.0;
    db.yh = 10.0;
    db.nodes = {
        {0, "a", 5.0, 5.0, 2.0, 2.0, false, false, "N"},
        {1, "b", 9.0, 5.0, 2.0, 2.0, false, false, "N"},
        {2, "target", 30.0, 5.0, 0.0, 0.0, true, true, "N"},
    };
    db.movable_ids = {0, 1};
    db.fixed_ids = {2};
    db.movable_area = 8.0;
    db.pins = {{0, 0.0, 0.0}, {1, 0.0, 0.0}, {2, 0.0, 0.0}};
    db.nets = {{0, "shared_target", 0, 3, 10.0}};
    db.node_pin_count = {1, 1, 1};
    db.node_pin_offsets = {0, 1, 2, 3};
    db.node_pin_indices = {0, 1, 2};
    ea::ExactOverlapDensity density(db, 4, 1, 1.0);
    const ea::DensityMetrics before_density =
        density.evaluate(0.0, 1.0, nullptr, nullptr);
    ea::RecoveryConfig config;
    config.sweeps = 1;
    config.line_search_steps = 1;
    config.step_bins = 0.5;
    config.hpwl_epsilon = 0.0;
    config.active_power = 1.0;
    config.degree_limit = 3;
    config.net_block = true;
    config.node_moves = false;
    config.net_block_degree_limit = 4;
    config.net_block_max_nodes = 2;
    config.net_block_max_blocks = 1;
    const double spacing = db.nodes[1].x - db.nodes[0].x;
    const ea::RecoveryStats stats =
        ea::recover_hpwl_under_overflow(db, density, config);
    require(stats.net_block_moves == 1 && stats.net_block_nodes == 2,
            "net block recovery must accept one rigid exact move");
    near(db.nodes[1].x - db.nodes[0].x, spacing, 1.0e-12,
         "net block move must preserve relative geometry");
    require(stats.final_hpwl < stats.initial_hpwl,
            "net block exact audit must reduce HPWL");
    require(stats.final_overflow <= before_density.overflow + 1.0e-12,
            "net block exact audit must not increase overlap");
}

void test_hierarchical_cluster_assignment() {
    ea::Database db;
    db.xl = 0.0;
    db.yl = 0.0;
    db.xh = 40.0;
    db.yh = 10.0;
    db.nodes = {
        {0, "a", 5.0, 4.0, 2.0, 2.0, false, false, "N"},
        {1, "b", 7.0, 6.0, 2.0, 2.0, false, false, "N"},
        {2, "source_macro", 5.0, 5.0, 10.0, 10.0, true, false, "N"},
        {3, "target", 35.0, 5.0, 0.0, 0.0, true, true, "N"},
    };
    db.movable_ids = {0, 1};
    db.fixed_ids = {2, 3};
    db.movable_area = 8.0;
    db.pins = {{0, 0.0, 0.0}, {1, 0.0, 0.0}, {3, 0.0, 0.0}};
    db.nets = {{0, "connected", 0, 3, 1.0}};
    db.node_pin_count = {1, 1, 0, 1};
    ea::ExactOverlapDensity density(db, 4, 1, 1.0);
    ea::DensityCoordinateConfig config;
    config.sweeps = 1;
    config.line_search_steps = 1;
    config.net_blocks = true;
    config.capacity_assignment = true;
    config.hierarchical_clusters = true;
    config.assignment_bins = 4;
    config.hierarchy_levels = 1;
    config.cluster_max_nodes = 4;
    config.cluster_candidate_bins = 4;
    config.cluster_max_clusters = 4;
    config.block_degree_limit = 4;
    config.hpwl_budget_fraction = 0.0;
    const double spacing_x = db.nodes[1].x - db.nodes[0].x;
    const double spacing_y = db.nodes[1].y - db.nodes[0].y;
    const double fixed_x = db.nodes[2].x;
    const ea::DensityCoordinateStats stats =
        ea::coordinate_descent_overlap(db, density, config);
    require(stats.moves >= 1,
            "hierarchical assignment must move one connected cluster");
    require(stats.final_overflow < stats.initial_overflow,
            "hierarchical assignment must reduce exact overlap");
    require(stats.final_hpwl < stats.initial_hpwl,
            "hierarchical assignment must use the topology-improving target");
    near(db.nodes[1].x - db.nodes[0].x, spacing_x, 1.0e-12,
         "hierarchical cluster must preserve relative x coordinates");
    near(db.nodes[1].y - db.nodes[0].y, spacing_y, 1.0e-12,
         "hierarchical cluster must preserve relative y coordinates");
    near(db.nodes[2].x, fixed_x, 1.0e-12,
         "hierarchical assignment must not move fixed macros");
}

}  // namespace

int main() {
    try {
        test_hpwl();
        test_fixed_macro_density();
        test_overlap_direction_and_exact_value();
        test_bin_field_prolongation();
        test_optimizer_factory();
        test_exact_batch_backtracking();
        test_exact_capacity_transport();
        test_recursive_bisection_fixed_capacity();
        test_constrained_hpwl_recovery();
        test_net_block_exact_recovery();
        test_hierarchical_cluster_assignment();
        test_density_neutral_swap_recovery();
        test_global_shape_permutation();
        test_net_disjoint_anchor_assignment();
        test_coarse_capacity_flow();
        std::cout << "all core tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
