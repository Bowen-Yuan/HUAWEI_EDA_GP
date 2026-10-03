#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <cmath>
#include <cfloat>
#include <algorithm>

using float_t = float;

struct Cell {
    int id;
    float_t x, y;         // current position
    float_t width, height;
    bool is_terminal;     // true = fixed IO pad
    float_t area;
};

struct Net {
    int id;
    std::vector<int> cell_ids;  // indices into cells array
    // Bookshelf pin offsets are relative to a cell centre, as in hpwl.pl.
    // Keeping them is required for parity with the supplied hpwl.pl script.
    std::vector<float_t> pin_offset_x;
    std::vector<float_t> pin_offset_y;
    float_t hpwl_x;
    float_t hpwl_y;
};

struct BinGrid {
    int nx, ny;
    float_t bin_w, bin_h;
    float_t chip_xl, chip_yl, chip_xh, chip_yh;
    float_t target_density;
    std::vector<float_t> density;     // nx * ny
    std::vector<float_t> overflow;    // nx * ny
    float_t total_overflow;
};

struct PlacementResult {
    float_t hpwl;
    float_t overflow;
    int iterations;
    double runtime_sec;
};

struct SolverConfig {
    int max_iters;           // Max outer iterations
    float_t step_scale;      // Base step size scaling
    float_t step_min;        // Minimum step size
    float_t step_max;        // Maximum step size
    float_t momentum;        // Heavy-ball momentum coefficient
    float_t lambda_init;     // Initial density penalty weight
    float_t lambda_growth;   // λ growth factor when overflow is high
    float_t lambda_decay;    // λ decay factor when overflow is low
    float_t target_overflow; // Target overflow ratio
    float_t density_target;  // Explicit ISPD density target
    float_t grad_clip;       // Gradient clipping threshold
    int p0_iters;            // Phase 0: HPWL-only iterations (default 150)
    int coarse_iters;        // Phase 1: iterations with coarse bins (default 100)
    int medium_iters;        // Phase 2: iterations with medium bins (default 150)
    float_t p1_lambda_init;  // Phase 1: initial λ for density penalty (default 0.005)
    bool use_bayesian_opt;   // Use Bayesian optimization for λ selection
};

struct CWTRConfig {
    int max_sweeps = 15;
    int candidates_per_dim = 5;
    float_t trust_radius_initial = 2.0f;
    float_t overflow_relax = 1.05f;
    float_t hpwl_improve_min = 0.9999f;
    int swap_attempts = 500;
};
