#include "types.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif

// Simple greedy row-based legalization
// Sorts cells by x-coordinate within each row and packs them left-to-right
// No smoothing — direct legalization of the non-smooth placement result

void legalize_placement(std::vector<Cell>& cells,
                         float_t /*chip_xl*/, float_t /*chip_yl*/,
                         float_t /*chip_xh*/, float_t /*chip_yh*/,
                         const std::string& scl_file)
{
    printf("[Legalizer] Starting legalization...\n");

    // Parse row definitions
    struct Row {
        float_t y_coord, height, sitewidth, sitespacing;
        float_t origin_x;
        int num_sites;
    };
    std::vector<Row> rows;

    {
        std::ifstream fin(scl_file);
        std::string line;
        Row cur_row;
        bool in_row = false;

        auto fast_split = [](const std::string& l, std::vector<std::string>& out) {
            out.clear();
            const char* s = l.c_str();
            while (*s) {
                while (*s == ' ' || *s == '\t' || *s == '\r') ++s;
                if (!*s) break;
                const char* start = s;
                while (*s && *s != ' ' && *s != '\t' && *s != '\r') ++s;
                out.emplace_back(start, s - start);
            }
        };

        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (line.find("CoreRow") != std::string::npos) {
                in_row = true;
                cur_row = Row{};
            } else if (line.find("End") != std::string::npos && in_row) {
                rows.push_back(cur_row);
                in_row = false;
            } else if (in_row) {
                if (line.find("Coordinate") != std::string::npos)
                    cur_row.y_coord = std::stof(line.substr(line.find(':') + 1));
                else if (line.find("Height") != std::string::npos)
                    cur_row.height = std::stof(line.substr(line.find(':') + 1));
                else if (line.find("Sitewidth") != std::string::npos)
                    cur_row.sitewidth = std::stof(line.substr(line.find(':') + 1));
                else if (line.find("Sitespacing") != std::string::npos)
                    cur_row.sitespacing = std::stof(line.substr(line.find(':') + 1));
                else if (line.find("SubrowOrigin") != std::string::npos) {
                    // Format: SubrowOrigin : <origin> NumSites : <numsites>
                    std::vector<std::string> tok;
                    fast_split(line, tok);
                    cur_row.origin_x = std::stof(tok[2]);
                    cur_row.num_sites = std::stoi(tok[5]);
                }
            }
        }
    }

    if (rows.empty()) {
        printf("[Legalizer] No rows found, skipping legalization\n");
        return;
    }

    // Assign each movable cell to the nearest row
    struct RowBin {
        std::vector<int> cell_ids;
        float_t cur_x;  // current packing position
    };
    std::vector<RowBin> row_bins(rows.size());

    for (size_t r = 0; r < rows.size(); ++r) {
        row_bins[r].cur_x = rows[r].origin_x;
    }

    // Collect movable cells for sorting
    std::vector<int> movable;
    for (int i = 0; i < (int)cells.size(); ++i)
        if (!cells[i].is_terminal)
            movable.push_back(i);

    // Sort cells by y-coordinate for row assignment
    std::sort(movable.begin(), movable.end(), [&](int a, int b) {
        return cells[a].y < cells[b].y;
    });

    // Assign cells to nearest rows
    #pragma omp parallel for schedule(dynamic, 256)
    for (int idx = 0; idx < (int)movable.size(); ++idx) {
        int i = movable[idx];
        float_t cy = cells[i].y;
        // Find nearest row by y-coordinate
        int best_r = 0;
        float_t best_dist = 1e30f;
        for (int r = 0; r < (int)rows.size(); ++r) {
            float_t ry = rows[r].y_coord + rows[r].height * 0.5f;
            float_t dist = std::abs(cy - ry);
            if (dist < best_dist) { best_dist = dist; best_r = r; }
        }
        // Snap to row y
        cells[i].y = rows[best_r].y_coord + rows[best_r].height * 0.5f;

        // Thread-safe addition to row bin
        #pragma omp critical
        {
            row_bins[best_r].cell_ids.push_back(i);
        }
    }

    // Pack cells within each row (left-to-right by x)
    for (size_t r = 0; r < rows.size(); ++r) {
        auto& bin = row_bins[r];
        if (bin.cell_ids.empty()) continue;

        // Sort cells within row by x-coordinate
        std::sort(bin.cell_ids.begin(), bin.cell_ids.end(), [&](int a, int b) {
            return cells[a].x < cells[b].x;
        });

        float_t row_x_end = rows[r].origin_x + rows[r].num_sites * rows[r].sitewidth;

        float_t cur_x = rows[r].origin_x;

        for (int cid : bin.cell_ids) {
            float_t cw = cells[cid].width;

            // Simple packing: place at cur_x, advance
            if (cur_x + cw > row_x_end) {
                // Row overflow — snap to start
                cur_x = rows[r].origin_x;
            }

            cells[cid].x = cur_x + cw * 0.5f;
            cur_x += cw;
        }
    }

    printf("[Legalizer] Legalized %zu cells into %zu rows\n",
           movable.size(), rows.size());
}
