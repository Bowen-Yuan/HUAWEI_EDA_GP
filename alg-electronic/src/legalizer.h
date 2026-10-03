#pragma once

#include "types.h"
#include <cstddef>
#include <string>
#include <vector>

struct LegalizationResult {
    bool success = false;
    std::size_t movable_cells = 0;
    std::size_t row_spans = 0;
    std::size_t free_segments = 0;
    double total_displacement = 0.0;
    double max_displacement = 0.0;
    double greedy_hpwl = 0.0;
    double abacus_hpwl = 0.0;
    double detailed_hpwl = 0.0;
    std::size_t detailed_swaps = 0;
    bool used_abacus = false;
    bool used_fallback = false;
    std::string method;
    std::string message;
};

enum class LegalizerMode {
    Legacy,
    DreamplaceGreedy,
    DreamplaceGreedyAbacus
};

struct LegalizerConfig {
    LegalizerMode mode = LegalizerMode::DreamplaceGreedyAbacus;
    double greedy_alpha = 0.5;
    int detailed_passes = 2;
    bool fallback_to_legacy = true;
};

LegalizationResult legalize_placement(std::vector<Cell>& cells,
                                      const std::vector<Net>& nets,
                                      float_t chip_xl, float_t chip_yl,
                                      float_t chip_xh, float_t chip_yh,
                                      const std::string& scl_file);

LegalizationResult legalize_placement(std::vector<Cell>& cells,
                                      const std::vector<Net>& nets,
                                      float_t chip_xl, float_t chip_yl,
                                      float_t chip_xh, float_t chip_yh,
                                      const std::string& scl_file,
                                      const LegalizerConfig& config);

LegalizationResult legalize_placement_legacy(std::vector<Cell>& cells,
                                             const std::vector<Net>& nets,
                                             float_t chip_xl, float_t chip_yl,
                                             float_t chip_xh, float_t chip_yh,
                                             const std::string& scl_file);
