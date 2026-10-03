#pragma once

#include "types.h"

void visualization_initialize(const VisualizationConfig& config,
                              const std::vector<Cell>& cells,
                              float_t chip_xl, float_t chip_yl,
                              float_t chip_xh, float_t chip_yh);
void visualization_snapshot(const VisualizationConfig& config, int iteration,
                            const char* phase, float_t hpwl, float_t overflow,
                            float_t density_penalty, float_t lambda,
                            const std::vector<Cell>& cells);
