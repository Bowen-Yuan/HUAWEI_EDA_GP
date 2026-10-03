#include "visualization.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace { bool initialized = false; }

void visualization_initialize(const VisualizationConfig& config,
                              const std::vector<Cell>& cells,
                              float_t chip_xl, float_t chip_yl,
                              float_t chip_xh, float_t chip_yh) {
    if (!config.enabled || initialized) return;
    std::filesystem::create_directories(config.output_dir);
    std::ofstream metadata(config.output_dir + "/metadata.txt");
    metadata << "nodes " << cells.size() << "\n"
             << "bounds " << chip_xl << " " << chip_yl << " " << chip_xh << " " << chip_yh << "\n";
    std::ofstream metrics(config.output_dir + "/metrics.csv");
    metrics << "iteration,phase,hpwl,overflow,density_penalty,lambda\n";
    initialized = true;
}

void visualization_snapshot(const VisualizationConfig& config, int iteration,
                            const char* phase, float_t hpwl, float_t overflow,
                            float_t density_penalty, float_t lambda,
                            const std::vector<Cell>& cells) {
    if (!config.enabled || config.every <= 0 || iteration % config.every != 0) return;
    char filename[64];
    std::snprintf(filename, sizeof(filename), "/snapshot_%06d.bin", iteration);
    std::ofstream snapshot(config.output_dir + filename, std::ios::binary);
    for (const Cell& cell : cells) {
        snapshot.write(reinterpret_cast<const char*>(&cell.x), sizeof(float_t));
        snapshot.write(reinterpret_cast<const char*>(&cell.y), sizeof(float_t));
    }
    std::ofstream metrics(config.output_dir + "/metrics.csv", std::ios::app);
    metrics << iteration << ',' << phase << ',' << std::setprecision(9) << hpwl << ','
            << overflow << ',' << density_penalty << ',' << lambda << '\n';
}
