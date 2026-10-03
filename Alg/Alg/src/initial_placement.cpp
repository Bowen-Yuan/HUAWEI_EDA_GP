#include "initial_placement.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <random>
#include <unordered_map>
#include <sstream>
#include <vector>
#include <cstring>

// Forward declaration
float_t compute_hpwl_only(const std::vector<Cell>&, const std::vector<Net>&);

// ---------------------------------------------------------------------------
// Parse .scl file for row information
// ---------------------------------------------------------------------------
std::vector<RowInfo> parse_scl_rows(const std::string& scl_file) {
    std::vector<RowInfo> rows;
    std::ifstream fin(scl_file);
    if (!fin) return rows;

    std::string line;
    RowInfo cur_row;
    bool in_row = false;

    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.find("CoreRow") != std::string::npos) {
            in_row = true;
            cur_row = RowInfo{};
        } else if (line.find("End") != std::string::npos && in_row) {
            rows.push_back(cur_row);
            in_row = false;
        } else if (in_row) {
            if (line.find("Coordinate") != std::string::npos)
                cur_row.coord = std::stof(line.substr(line.find(':') + 1));
            else if (line.find("Height") != std::string::npos)
                cur_row.height = std::stof(line.substr(line.find(':') + 1));
            else if (line.find("Sitewidth") != std::string::npos)
                cur_row.sitewidth = std::stof(line.substr(line.find(':') + 1));
            else if (line.find("Sitespacing") != std::string::npos)
                cur_row.sitespacing = std::stof(line.substr(line.find(':') + 1));
            else if (line.find("SubrowOrigin") != std::string::npos) {
                // Format: SubrowOrigin : <origin> NumSites : <numsites>
                const char* s = line.c_str();
                // skip "SubrowOrigin : "
                while (*s && *s != ':') ++s;
                if (*s) ++s; // skip ':'
                while (*s == ' ') ++s;
                cur_row.origin = std::stof(std::string(s, strchr(s, ' ') - s));
                // find "NumSites :"
                const char* ns = strstr(s, "NumSites");
                if (ns) {
                    ns += 8; // skip "NumSites"
                    while (*ns == ' ' || *ns == ':') ++ns;
                    cur_row.numsites = std::stoi(ns);
                }
            }
        }
    }
    return rows;
}

// ---------------------------------------------------------------------------
// Helper: collect movable cell IDs
// ---------------------------------------------------------------------------
static std::vector<int> get_movable_ids(const std::vector<Cell>& cells) {
    std::vector<int> ids;
    for (int i = 0; i < (int)cells.size(); ++i)
        if (!cells[i].is_terminal)
            ids.push_back(i);
    return ids;
}

// ---------------------------------------------------------------------------
// Strategy implementations
// ---------------------------------------------------------------------------

// ---- Random placement ----
static void init_random(std::vector<Cell>& cells, float_t xl, float_t yl,
                         float_t xh, float_t yh) {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<float_t> dist_x(xl, xh);
    std::uniform_real_distribution<float_t> dist_y(yl, yh);

    for (auto& c : cells) {
        if (c.is_terminal) continue;
        c.x = dist_x(rng);
        c.y = dist_y(rng);
    }
    printf("  [Random] Cells uniformly distributed in chip area\n");
}

// ---- Uniform grid ----
static void init_uniform_grid(std::vector<Cell>& cells, float_t xl, float_t yl,
                               float_t xh, float_t yh) {
    auto movable = get_movable_ids(cells);
    int Nm = (int)movable.size();
    float_t cw = xh - xl, ch = yh - yl;

    // Sort by width for better packing
    std::sort(movable.begin(), movable.end(), [&](int a, int b) {
        return cells[a].width > cells[b].width;
    });

    int cols = (int)std::sqrt((double)Nm * cw / ch);
    cols = std::max(1, cols);
    float_t sp_x = cw / cols;
    float_t sp_y = ch / ((Nm + cols - 1) / cols);

    for (int idx = 0; idx < Nm; ++idx) {
        int i = movable[idx];
        int row = idx / cols, col = idx % cols;
        cells[i].x = xl + (col + 0.5f) * sp_x;
        cells[i].y = yl + (row + 0.5f) * sp_y;
        cells[i].x = std::max(xl + cells[i].width*0.5f,
                               std::min(xh - cells[i].width*0.5f, cells[i].x));
        cells[i].y = std::max(yl + cells[i].height*0.5f,
                               std::min(yh - cells[i].height*0.5f, cells[i].y));
    }
    printf("  [Uniform Grid] %d cells in %d cols\n", Nm, cols);
}

// ---- Row-ordered placement ----
static void init_row_ordered(std::vector<Cell>& cells, float_t xl, float_t yl,
                              float_t xh, float_t yh,
                              const std::vector<RowInfo>& rows) {
    if (rows.empty()) {
        printf("  [Row-Ordered] No row info, falling back to uniform grid\n");
        init_uniform_grid(cells, xl, yl, xh, yh);
        return;
    }

    auto movable = get_movable_ids(cells);

    // Sort cells by height for row matching
    std::sort(movable.begin(), movable.end(), [&](int a, int b) {
        return cells[a].height > cells[b].height;
    });

    // Assign cells to rows round-robin
    int n_rows = (int)rows.size();
    float_t total_width = 0;
    for (auto& r : rows) total_width += r.numsites * r.sitewidth;

    size_t cell_idx = 0;
    for (int ri = 0; ri < n_rows && cell_idx < movable.size(); ++ri) {
        const auto& row = rows[ri];
        float_t row_x_start = row.origin;
        float_t row_width = row.numsites * row.sitewidth;
        float_t row_y = row.coord + row.height * 0.5f;

        int cells_in_row = (int)((float)movable.size() * row_width / total_width);
        cells_in_row = std::max(1, cells_in_row);

        for (int c = 0; c < cells_in_row && cell_idx < movable.size(); ++c, ++cell_idx) {
            int id = movable[cell_idx];
            float_t frac = (c + 0.5f) / cells_in_row;
            cells[id].x = row_x_start + frac * row_width;
            cells[id].y = row_y;
            cells[id].x = std::max(xl + cells[id].width*0.5f,
                                    std::min(xh - cells[id].width*0.5f, cells[id].x));
        }
    }

    printf("  [Row-Ordered] %zu cells placed across %d rows\n", cell_idx, n_rows);
}

// ---- Scaled ePlace: take ePlace coords and expand from center ----
static void init_scaled_eplace(std::vector<Cell>& cells, float_t xl, float_t yl,
                                float_t xh, float_t yh, float_t scale_factor) {
    float_t cx = (xl + xh) * 0.5f;
    float_t cy = (yl + yh) * 0.5f;
    int count = 0;

    for (auto& c : cells) {
        if (c.is_terminal) continue;
        if (c.x == 0 && c.y == 0) continue; // was never placed
        float_t dx = c.x - cx;
        float_t dy = c.y - cy;
        c.x = cx + dx * scale_factor;
        c.y = cy + dy * scale_factor;
        // Allow cells to go slightly outside chip (will be clamped by solver)
        c.x = std::max(xl - c.width, std::min(xh + c.width, c.x));
        c.y = std::max(yl - c.height, std::min(yh + c.height, c.y));
        ++count;
    }
    printf("  [Scaled ePlace] %d cells expanded by factor %.2f from center\n",
           count, scale_factor);
}

// ---- Gaussian around center ----
static void init_gaussian_center(std::vector<Cell>& cells, float_t xl, float_t yl,
                                  float_t xh, float_t yh) {
    std::mt19937 rng(67890);
    float_t cx = (xl + xh) * 0.5f;
    float_t cy = (yl + yh) * 0.5f;
    float_t sx = (xh - xl) * 0.25f;
    float_t sy = (yh - yl) * 0.25f;
    std::normal_distribution<float_t> dist_x(0.0f, sx);
    std::normal_distribution<float_t> dist_y(0.0f, sy);

    for (auto& c : cells) {
        if (c.is_terminal) continue;
        c.x = cx + dist_x(rng);
        c.y = cy + dist_y(rng);
        c.x = std::max(xl, std::min(xh, c.x));
        c.y = std::max(yl, std::min(yh, c.y));
    }
    printf("  [Gaussian Center] Cells around (%.0f,%.0f) with sigma=(%.0f,%.0f)\n",
           cx, cy, sx, sy);
}

// ---- Quadrant bins: use net connectivity to assign cells ----
static void init_quadrant_bins(std::vector<Cell>& cells,
                                const std::vector<Net>& nets,
                                float_t xl, float_t yl, float_t xh, float_t yh) {
    auto movable = get_movable_ids(cells);

    // Count terminal connectivity per cell
    std::vector<int> term_conn(cells.size(), 0);
    std::vector<float_t> term_avg_x(cells.size(), 0.0f);
    std::vector<float_t> term_avg_y(cells.size(), 0.0f);

    for (const auto& net : nets) {
        float_t tx_sum = 0, ty_sum = 0;
        int tcount = 0;
        for (int cid : net.cell_ids) {
            if (cells[cid].is_terminal) {
                tx_sum += cells[cid].x;
                ty_sum += cells[cid].y;
                ++tcount;
            }
        }
        if (tcount > 0) {
            float_t tax = tx_sum / tcount;
            float_t tay = ty_sum / tcount;
            for (int cid : net.cell_ids) {
                if (!cells[cid].is_terminal) {
                    term_conn[cid]++;
                    term_avg_x[cid] += tax;
                    term_avg_y[cid] += tay;
                }
            }
        }
    }

    // For cells with terminal connections, place near terminal centroid
    std::mt19937 rng(11111);
    std::uniform_real_distribution<float_t> jitter(-500, 500);

    int placed_by_term = 0, placed_random = 0;
    for (auto& c : cells) {
        if (c.is_terminal) continue;
        int id = c.id;
        if (term_conn[id] > 0) {
            c.x = term_avg_x[id] / term_conn[id] + jitter(rng);
            c.y = term_avg_y[id] / term_conn[id] + jitter(rng);
            c.x = std::max(xl, std::min(xh, c.x));
            c.y = std::max(yl, std::min(yh, c.y));
            ++placed_by_term;
        } else {
            // fallback: random
            std::uniform_real_distribution<float_t> rx(xl, xh);
            std::uniform_real_distribution<float_t> ry(yl, yh);
            c.x = rx(rng);
            c.y = ry(rng);
            ++placed_random;
        }
    }
    printf("  [Quadrant Bins] %d cells by terminal centroid, %d random fallback\n",
           placed_by_term, placed_random);
}

// ---------------------------------------------------------------------------
// Main dispatch function
// ---------------------------------------------------------------------------
std::string apply_initial_placement(
    InitStrategy strategy,
    std::vector<Cell>& cells,
    const std::vector<Net>& nets,
    float_t chip_xl, float_t chip_yl,
    float_t chip_xh, float_t chip_yh,
    const std::string& bench_name)
{
    switch (strategy) {
    case InitStrategy::EPLACE_QP:
        // Already loaded by parser — no-op
        printf("[Init] Using ePlace QP placement (loaded from file)\n");
        return "ePlace_QP";

    case InitStrategy::UNIFORM_GRID:
        init_uniform_grid(cells, chip_xl, chip_yl, chip_xh, chip_yh);
        return "Uniform_Grid";

    case InitStrategy::RANDOM:
        init_random(cells, chip_xl, chip_yl, chip_xh, chip_yh);
        return "Random";

    case InitStrategy::ROW_ORDERED: {
        std::string scl_file = bench_name + ".scl";
        auto rows = parse_scl_rows(scl_file);
        init_row_ordered(cells, chip_xl, chip_yl, chip_xh, chip_yh, rows);
        return "Row_Ordered";
    }

    case InitStrategy::SCALED_EPLACE:
        // Expand ePlace coords by 1.3x from center to reduce overlap
        init_scaled_eplace(cells, chip_xl, chip_yl, chip_xh, chip_yh, 1.30f);
        return "Scaled_ePlace_x1.3";

    case InitStrategy::GAUSSIAN_CENTER:
        init_gaussian_center(cells, chip_xl, chip_yl, chip_xh, chip_yh);
        return "Gaussian_Center";

    case InitStrategy::QUADRANT_BINS:
        init_quadrant_bins(cells, nets, chip_xl, chip_yl, chip_xh, chip_yh);
        return "Quadrant_Bins";

    default:
        printf("[Init] Unknown strategy, keeping current positions\n");
        return "Unknown";
    }
}
