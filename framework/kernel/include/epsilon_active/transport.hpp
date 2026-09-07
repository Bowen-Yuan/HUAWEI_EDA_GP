#pragma once

#include "epsilon_active/density.hpp"
#include "epsilon_active/types.hpp"

#include <string>

namespace ea {

enum class TransportDestinationMode { Nearest, Auction, Hilbert };

TransportDestinationMode parse_transport_destination_mode(const std::string& name);
const char* transport_destination_mode_name(TransportDestinationMode mode) noexcept;

struct TransportConfig {
    int rounds = 0;
    int max_moves = 250000;
    int max_source_bins = 64;
    int candidate_lookahead = 8;
    int group_size = 1;
    int group_degree_limit = 32;
    int group_destination_radius = 0;
    bool group_collective_target = false;
    bool weighted_grouping = false;
    bool preserve_bin_offset = false;
    bool rigid_group_displacement = false;
    int rigid_rounds = -1;
    bool atomic_groups = false;
    bool atomic_source_active = false;
    int atomic_rounds = -1;
    bool atomic_componentwise_capacity = false;
    int atomic_capacity_candidates = 1;
    bool atomic_source_batch = false;
    Real atomic_min_gain_ratio = 0.0;
    bool atomic_fallback_individual = false;
    int identity_exchange_passes = 0;
    int identity_exchange_candidates = 8;
    int identity_exchange_radius = 16;
    TransportDestinationMode destination_mode = TransportDestinationMode::Nearest;
    int auction_candidates = 128;
    int auction_shortlist = 8;
    int auction_hpwl_degree_limit = 100;
    Real auction_density_price_weight = 1.0;
    int auction_rounds = -1;
    bool connectivity_order = false;
    int connectivity_degree_limit = 32;
    int hilbert_rounds = 1;
    int hilbert_window = 128;
    // Optional cumulative guard per transport round. Negative disables it;
    // nonnegative values are fractional HPWL budgets above round start.
    Real global_hpwl_budget = -1.0;
};

struct TransportStats {
    int rounds = 0;
    int moves = 0;
    int groups = 0;
    int atomic_attempts = 0;
    int atomic_groups = 0;
    int atomic_fallbacks = 0;
    int atomic_capacity_rejections = 0;
    int atomic_bid_trials = 0;
    int atomic_feasible_bids = 0;
    int atomic_rescued_groups = 0;
    int atomic_batch_sources = 0;
    int atomic_batch_groups = 0;
    int atomic_batch_bids = 0;
    int atomic_batch_rejections = 0;
    int atomic_batch_commits = 0;
    int identity_exchange_attempts = 0;
    int identity_exchange_accepted = 0;
    int identity_exchange_permuted_nodes = 0;
    Real identity_initial_hpwl = 0.0;
    Real identity_final_hpwl = 0.0;
    Real identity_initial_overflow = 0.0;
    Real identity_final_overflow = 0.0;
    double identity_wall_seconds = 0.0;
    Real initial_hpwl = 0.0;
    Real initial_overflow = 0.0;
    Real final_hpwl = 0.0;
    Real final_overflow = 0.0;
    double wall_seconds = 0.0;
};

TransportStats transport_excess_to_capacity(
    Database& db, ExactOverlapDensity& density, const TransportConfig& config);

}  // namespace ea
