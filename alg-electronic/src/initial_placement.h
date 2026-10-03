#pragma once
#include "types.h"
#include <string>

// ============================================================================
// Initial Placement Strategies
//
// Different strategies for setting the starting positions of movable cells
// before the main optimization begins. The choice of initial placement can
// dramatically affect convergence speed and final solution quality.
// ============================================================================

enum class InitStrategy {
    EPLACE_QP,       // Load from ePlace .pl file (default) — best HPWL, high overlap
    UNIFORM_GRID,    // Uniform sqrt(N) grid — low overflow, terrible HPWL
    RANDOM,          // Uniform random within chip area
    ROW_ORDERED,     // Pack cells in row order from .scl — respects row structure
    SCALED_EPLACE,   // Take ePlace coords and expand from center (reduce overlap)
    GAUSSIAN_CENTER, // Place around chip center with Gaussian noise
    DREAMPLACE_CENTER, // Center with 0.1% Gaussian noise, as in DREAMPlace
    QUADRANT_BINS,   // Assign cells to quadrants based on net connectivity
};

struct RowInfo {
    float_t coord;      // y-coordinate of row
    float_t height;     // row height
    float_t sitewidth;  // site width
    float_t sitespacing;// site spacing
    float_t origin;     // x-origin of subrow
    int numsites;       // number of sites in subrow
};

// Parse .scl file to extract row layout information
std::vector<RowInfo> parse_scl_rows(const std::string& scl_file);

// Apply the selected initial placement strategy
// Returns the name of the strategy applied (for logging)
std::string apply_initial_placement(
    InitStrategy strategy,
    std::vector<Cell>& cells,
    const std::vector<Net>& nets,
    float_t chip_xl, float_t chip_yl,
    float_t chip_xh, float_t chip_yh,
    const std::string& bench_name);
