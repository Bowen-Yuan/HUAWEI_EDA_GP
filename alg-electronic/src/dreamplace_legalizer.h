#pragma once

#include "legalizer.h"

LegalizationResult dreamplace_legalize_placement(
    std::vector<Cell>& cells, const std::vector<Net>& nets,
    float_t chip_xl, float_t chip_yl, float_t chip_xh, float_t chip_yh,
    const std::string& scl_file, const LegalizerConfig& config);
