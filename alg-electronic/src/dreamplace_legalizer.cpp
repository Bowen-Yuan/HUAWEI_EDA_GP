// Adaptation of DREAMPlace's CPU greedy and Abacus legalization flow.
// DREAMPlace is distributed under the BSD 3-Clause License; see
// THIRD_PARTY_NOTICES.md. Coordinates are converted from this placer's cell
// centers to DREAMPlace-style lower-left locations inside this file.

#include "dreamplace_legalizer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

float_t compute_hpwl_only(const std::vector<Cell>& cells,
                          const std::vector<Net>& nets);

namespace {

constexpr double kEps = 1.0e-5;

struct Row {
    double y = 0.0;
    double height = 0.0;
    double origin = 0.0;
    double xh = 0.0;
    double site = 1.0;
};

struct Blank {
    double xl = 0.0;
    double xh = 0.0;
};

struct ParsedRow {
    double y = 0.0;
    double height = 0.0;
    double site_width = 0.0;
    double site_spacing = 0.0;
    std::vector<std::pair<double, int>> subrows;
};

std::vector<std::string> split_ws(const std::string& line) {
    std::vector<std::string> out;
    const char* p = line.c_str();
    while (*p) {
        while (*p == ' ' || *p == '\t' || *p == '\r') ++p;
        if (!*p) break;
        const char* begin = p;
        while (*p && *p != ' ' && *p != '\t' && *p != '\r') ++p;
        out.emplace_back(begin, p - begin);
    }
    return out;
}

bool parse_rows(const std::string& scl_file, std::vector<Row>& rows,
                std::string& error) {
    std::ifstream fin(scl_file);
    if (!fin) {
        error = "cannot open SCL file: " + scl_file;
        return false;
    }
    ParsedRow parsed;
    bool in_row = false;
    std::string line;
    while (std::getline(fin, line)) {
        const auto tok = split_ws(line);
        if (tok.empty() || tok[0] == "#" || tok[0] == "UCLA") continue;
        if (tok[0] == "CoreRow") {
            parsed = ParsedRow{};
            in_row = true;
            continue;
        }
        if (!in_row) continue;
        if (tok[0] == "End") {
            const double site = parsed.site_spacing > 0.0
                ? parsed.site_spacing : parsed.site_width;
            if (parsed.height <= 0.0 || site <= 0.0 || parsed.subrows.empty()) {
                error = "incomplete CoreRow in SCL file";
                return false;
            }
            for (const auto& subrow : parsed.subrows) {
                if (subrow.second > 0) {
                    rows.push_back(Row{parsed.y, parsed.height, subrow.first,
                                       subrow.first + subrow.second * site, site});
                }
            }
            in_row = false;
        } else if (tok[0] == "Coordinate" && tok.size() >= 3) {
            parsed.y = std::stod(tok[2]);
        } else if (tok[0] == "Height" && tok.size() >= 3) {
            parsed.height = std::stod(tok[2]);
        } else if (tok[0] == "Sitewidth" && tok.size() >= 3) {
            parsed.site_width = std::stod(tok[2]);
        } else if (tok[0] == "Sitespacing" && tok.size() >= 3) {
            parsed.site_spacing = std::stod(tok[2]);
        } else if (tok[0] == "SubrowOrigin" && tok.size() >= 6) {
            parsed.subrows.emplace_back(std::stod(tok[2]), std::stoi(tok[5]));
        }
    }
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
        return a.y < b.y || (a.y == b.y && a.origin < b.origin);
    });
    if (rows.empty()) error = "no placement rows found in SCL file";
    return !rows.empty();
}

double snap_up(double x, double origin, double site) {
    return origin + std::ceil((x - origin) / site - kEps) * site;
}

double snap_down(double x, double origin, double site) {
    return origin + std::floor((x - origin) / site + kEps) * site;
}

std::vector<std::vector<Blank>> build_blanks(const std::vector<Cell>& cells,
                                              const std::vector<Row>& rows) {
    std::vector<std::vector<Blank>> all(rows.size());
    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
        const Row& row = rows[r];
        std::vector<Blank> blocked;
        for (const Cell& cell : cells) {
            if (!cell.is_terminal) continue;
            const double yl = cell.y - 0.5 * cell.height;
            const double yh = cell.y + 0.5 * cell.height;
            if (yl >= row.y + row.height - kEps || yh <= row.y + kEps) continue;
            const double xl = std::max(row.origin,
                                       static_cast<double>(cell.x - 0.5f * cell.width));
            const double xh = std::min(row.xh,
                                       static_cast<double>(cell.x + 0.5f * cell.width));
            if (xh > xl + kEps) blocked.push_back({xl, xh});
        }
        std::sort(blocked.begin(), blocked.end(), [](const Blank& a, const Blank& b) {
            return a.xl < b.xl || (a.xl == b.xl && a.xh < b.xh);
        });
        std::vector<Blank> merged;
        for (const Blank& b : blocked) {
            if (merged.empty() || b.xl > merged.back().xh + kEps) {
                merged.push_back(b);
            } else {
                merged.back().xh = std::max(merged.back().xh, b.xh);
            }
        }
        double cursor = row.origin;
        auto append = [&](double lo, double hi) {
            lo = snap_up(lo, row.origin, row.site);
            hi = snap_down(hi, row.origin, row.site);
            if (hi > lo + kEps) all[r].push_back({lo, hi});
        };
        for (const Blank& b : merged) {
            append(cursor, b.xl);
            cursor = std::max(cursor, b.xh);
        }
        append(cursor, row.xh);
    }
    return all;
}

bool greedy_sweep(std::vector<Cell>& cells, const std::vector<Cell>& initial,
                  const std::vector<Row>& rows,
                  const std::vector<std::vector<Blank>>& initial_blanks,
                  double alpha, bool left_to_right, std::string& error) {
    std::vector<std::vector<Blank>> blanks = initial_blanks;
    std::vector<int> movable;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        if (!cells[i].is_terminal) movable.push_back(i);
    }
    std::sort(movable.begin(), movable.end(), [&](int i, int j) {
        const double ci = initial[i].x;
        const double cj = initial[j].x;
        if (ci != cj) return left_to_right ? ci < cj : ci > cj;
        if (initial[i].y != initial[j].y) {
            return left_to_right ? initial[i].y < initial[j].y
                                 : initial[i].y > initial[j].y;
        }
        return i < j;
    });

    for (int id : movable) {
        const double width = cells[id].width;
        const double height = cells[id].height;
        const double target_center_x = alpha * initial[id].x +
                                       (1.0 - alpha) * cells[id].x;
        const double target_center_y = alpha * initial[id].y +
                                       (1.0 - alpha) * cells[id].y;
        int best_row = -1;
        int best_blank = -1;
        double best_left = 0.0;
        double best_cost = std::numeric_limits<double>::infinity();

        for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
            const Row& row = rows[r];
            if (height > row.height + kEps) continue;
            const double placed_y = row.y + 0.5 * height;
            const double y_cost = std::abs(placed_y - target_center_y);
            if (y_cost > best_cost + row.height) continue;
            for (int b = 0; b < static_cast<int>(blanks[r].size()); ++b) {
                const Blank& blank = blanks[r][b];
                const double rounded_width = snap_up(width, 0.0, row.site);
                const double available = blank.xh - blank.xl;
                if (available + kEps < rounded_width) continue;
                double left = snap_down(target_center_x - 0.5 * width,
                                        row.origin, row.site);
                left = std::clamp(left, blank.xl, blank.xh - rounded_width);
                const double tolerance = std::min(4.0 * rounded_width,
                                                  available / 4.0);
                if (left <= blank.xl + tolerance) left = blank.xl;
                else if (left + rounded_width >= blank.xh - tolerance) {
                    left = blank.xh - rounded_width;
                }
                const double cost = y_cost +
                    std::abs((left + 0.5 * width) - target_center_x);
                if (cost < best_cost - kEps) {
                    best_row = r;
                    best_blank = b;
                    best_left = left;
                    best_cost = cost;
                }
            }
        }
        if (best_row < 0) {
            error = "DREAMPlace greedy cannot place movable cell o" +
                    std::to_string(cells[id].id);
            return false;
        }

        const Row& row = rows[best_row];
        const double rounded_width = snap_up(width, 0.0, row.site);
        std::vector<Blank>& row_blanks = blanks[best_row];
        const Blank selected = row_blanks[best_blank];
        const double right = best_left + rounded_width;
        if (best_left <= selected.xl + kEps && right >= selected.xh - kEps) {
            row_blanks.erase(row_blanks.begin() + best_blank);
        } else if (best_left <= selected.xl + kEps) {
            row_blanks[best_blank].xl = right;
        } else if (right >= selected.xh - kEps) {
            row_blanks[best_blank].xh = best_left;
        } else {
            row_blanks[best_blank].xh = best_left;
            row_blanks.insert(row_blanks.begin() + best_blank + 1,
                              Blank{right, selected.xh});
        }
        cells[id].x = static_cast<float_t>(best_left + 0.5 * width);
        cells[id].y = static_cast<float_t>(row.y + 0.5 * height);
    }
    return true;
}

bool dreamplace_greedy(std::vector<Cell>& cells,
                       const std::vector<Row>& rows,
                       const std::vector<std::vector<Blank>>& blanks,
                       double alpha, std::string& error) {
    const std::vector<Cell> initial = cells;
    // DREAMPlace performs a right-to-left sweep followed by a left-to-right
    // sweep. Blanks are rebuilt between sweeps while positions are retained.
    return greedy_sweep(cells, initial, rows, blanks, alpha, false, error) &&
           greedy_sweep(cells, initial, rows, blanks, alpha, true, error);
}

struct Cluster {
    int first = 0;
    int last = 0;
    double e = 0.0;
    double q = 0.0;
    double width = 0.0;
    double x = 0.0;
};

void abacus_segment(std::vector<Cell>& cells, const std::vector<Cell>& initial,
                    const Row& row, const Blank& segment,
                    std::vector<int> ids) {
    if (ids.empty()) return;
    std::sort(ids.begin(), ids.end(), [&](int a, int b) {
        const double ax = cells[a].x - 0.5 * cells[a].width;
        const double bx = cells[b].x - 0.5 * cells[b].width;
        return ax < bx || (ax == bx && a < b);
    });
    std::vector<double> widths(ids.size());
    std::vector<Cluster> clusters;
    clusters.reserve(ids.size());
    for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
        widths[i] = snap_up(cells[ids[i]].width, 0.0, row.site);
        const double target = initial[ids[i]].x - 0.5 * initial[ids[i]].width;
        Cluster c{i, i, 1.0, target, widths[i], 0.0};
        c.x = std::clamp(c.q / c.e, segment.xl, segment.xh - c.width);
        clusters.push_back(c);
        while (clusters.size() >= 2) {
            Cluster& right = clusters.back();
            Cluster& left = clusters[clusters.size() - 2];
            if (left.x + left.width <= right.x + kEps) break;
            left.last = right.last;
            left.q += right.q - right.e * left.width;
            left.e += right.e;
            left.width += right.width;
            left.x = std::clamp(left.q / left.e, segment.xl,
                                segment.xh - left.width);
            clusters.pop_back();
        }
    }
    for (const Cluster& c : clusters) {
        double left = snap_down(c.x, row.origin, row.site);
        left = std::clamp(left, segment.xl, segment.xh - c.width);
        for (int i = c.first; i <= c.last; ++i) {
            const int id = ids[i];
            cells[id].x = static_cast<float_t>(left + 0.5 * cells[id].width);
            cells[id].y = static_cast<float_t>(row.y + 0.5 * cells[id].height);
            left += widths[i];
        }
    }
}

bool run_abacus(std::vector<Cell>& cells, const std::vector<Cell>& initial,
                const std::vector<Row>& rows,
                const std::vector<std::vector<Blank>>& segments,
                std::string& error) {
    std::vector<std::vector<std::vector<int>>> assignments(rows.size());
    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
        assignments[r].resize(segments[r].size());
    }
    for (int id = 0; id < static_cast<int>(cells.size()); ++id) {
        if (cells[id].is_terminal) continue;
        const double left = cells[id].x - 0.5 * cells[id].width;
        const double bottom = cells[id].y - 0.5 * cells[id].height;
        bool found = false;
        for (int r = 0; r < static_cast<int>(rows.size()) && !found; ++r) {
            if (std::abs(bottom - rows[r].y) > kEps) continue;
            for (int s = 0; s < static_cast<int>(segments[r].size()); ++s) {
                if (left >= segments[r][s].xl - kEps &&
                    left + cells[id].width <= segments[r][s].xh + kEps) {
                    assignments[r][s].push_back(id);
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            error = "cannot map greedy cell o" + std::to_string(cells[id].id) +
                    " to an Abacus row segment";
            return false;
        }
    }
    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
        for (int s = 0; s < static_cast<int>(segments[r].size()); ++s) {
            abacus_segment(cells, initial, rows[r], segments[r][s],
                           std::move(assignments[r][s]));
        }
    }
    return true;
}

bool check_legality(const std::vector<Cell>& cells,
                    const std::vector<Row>& rows,
                    const std::vector<std::vector<Blank>>& segments,
                    std::string& error) {
    std::vector<std::vector<std::pair<double, double>>> occupied(rows.size());
    for (int id = 0; id < static_cast<int>(cells.size()); ++id) {
        if (cells[id].is_terminal) continue;
        const double left = cells[id].x - 0.5 * cells[id].width;
        const double right = cells[id].x + 0.5 * cells[id].width;
        const double bottom = cells[id].y - 0.5 * cells[id].height;
        bool found = false;
        for (int r = 0; r < static_cast<int>(rows.size()) && !found; ++r) {
            if (std::abs(bottom - rows[r].y) > kEps ||
                cells[id].height > rows[r].height + kEps) continue;
            const double snapped = snap_down(left, rows[r].origin, rows[r].site);
            if (std::abs(snapped - left) > kEps) continue;
            for (const Blank& segment : segments[r]) {
                if (left >= segment.xl - kEps && right <= segment.xh + kEps) {
                    occupied[r].push_back({left, right});
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            error = "cell o" + std::to_string(cells[id].id) +
                    " is not row/site/obstacle legal";
            return false;
        }
    }
    for (auto& row : occupied) {
        std::sort(row.begin(), row.end());
        for (int i = 1; i < static_cast<int>(row.size()); ++i) {
            if (row[i].first < row[i - 1].second - kEps) {
                error = "movable-cell overlap remains after legalization";
                return false;
            }
        }
    }
    return true;
}

void displacement(const std::vector<Cell>& cells,
                  const std::vector<Cell>& initial,
                  LegalizationResult& result) {
    result.total_displacement = 0.0;
    result.max_displacement = 0.0;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        if (cells[i].is_terminal) continue;
        const double d = std::hypot(static_cast<double>(cells[i].x - initial[i].x),
                                    static_cast<double>(cells[i].y - initial[i].y));
        result.total_displacement += d;
        result.max_displacement = std::max(result.max_displacement, d);
    }
}

double net_hpwl(const std::vector<Cell>& cells, const Net& net) {
    if (net.cell_ids.size() < 2) return 0.0;
    double min_x = cells[net.cell_ids[0]].x;
    double max_x = min_x;
    double min_y = cells[net.cell_ids[0]].y;
    double max_y = min_y;
    for (std::size_t i = 1; i < net.cell_ids.size(); ++i) {
        const Cell& cell = cells[net.cell_ids[i]];
        min_x = std::min(min_x, static_cast<double>(cell.x));
        max_x = std::max(max_x, static_cast<double>(cell.x));
        min_y = std::min(min_y, static_cast<double>(cell.y));
        max_y = std::max(max_y, static_cast<double>(cell.y));
    }
    return max_x - min_x + max_y - min_y;
}

std::size_t detailed_adjacent_reorder(
    std::vector<Cell>& cells, const std::vector<Net>& nets,
    const std::vector<Row>& rows,
    const std::vector<std::vector<Blank>>& segments, int passes) {
    if (passes <= 0 || nets.empty()) return 0;
    std::vector<std::vector<int>> cell_nets(cells.size());
    for (int n = 0; n < static_cast<int>(nets.size()); ++n) {
        for (int id : nets[n].cell_ids) cell_nets[id].push_back(n);
    }
    std::vector<std::vector<std::vector<int>>> assignments(rows.size());
    for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
        assignments[r].resize(segments[r].size());
    }
    for (int id = 0; id < static_cast<int>(cells.size()); ++id) {
        if (cells[id].is_terminal) continue;
        const double left = cells[id].x - 0.5 * cells[id].width;
        const double bottom = cells[id].y - 0.5 * cells[id].height;
        for (int r = 0; r < static_cast<int>(rows.size()); ++r) {
            if (std::abs(bottom - rows[r].y) > kEps) continue;
            for (int s = 0; s < static_cast<int>(segments[r].size()); ++s) {
                if (left >= segments[r][s].xl - kEps &&
                    left + cells[id].width <= segments[r][s].xh + kEps) {
                    assignments[r][s].push_back(id);
                    break;
                }
            }
            break;
        }
    }
    for (auto& row : assignments) {
        for (auto& ids : row) {
            std::sort(ids.begin(), ids.end(), [&](int a, int b) {
                return cells[a].x < cells[b].x ||
                       (cells[a].x == cells[b].x && a < b);
            });
        }
    }

    std::vector<int> marks(nets.size(), 0);
    int token = 0;
    std::vector<int> affected;
    std::size_t accepted = 0;
    for (int pass = 0; pass < passes; ++pass) {
        std::size_t pass_accepted = 0;
        for (auto& row : assignments) {
            for (auto& ids : row) {
                const int parity = pass & 1;
                for (int i = parity; i + 1 < static_cast<int>(ids.size()); i += 2) {
                    const int a = ids[i];
                    const int b = ids[i + 1];
                    const double a_left = cells[a].x - 0.5 * cells[a].width;
                    const double b_left = cells[b].x - 0.5 * cells[b].width;
                    const double gap = b_left - (a_left + cells[a].width);
                    if (gap < -kEps) continue;

                    if (++token == std::numeric_limits<int>::max()) {
                        std::fill(marks.begin(), marks.end(), 0);
                        token = 1;
                    }
                    affected.clear();
                    for (int n : cell_nets[a]) {
                        if (marks[n] != token) {
                            marks[n] = token;
                            affected.push_back(n);
                        }
                    }
                    for (int n : cell_nets[b]) {
                        if (marks[n] != token) {
                            marks[n] = token;
                            affected.push_back(n);
                        }
                    }
                    double before = 0.0;
                    for (int n : affected) before += net_hpwl(cells, nets[n]);

                    const float_t old_ax = cells[a].x;
                    const float_t old_bx = cells[b].x;
                    cells[b].x = static_cast<float_t>(
                        a_left + 0.5 * cells[b].width);
                    cells[a].x = static_cast<float_t>(
                        a_left + cells[b].width + gap + 0.5 * cells[a].width);
                    double after = 0.0;
                    for (int n : affected) after += net_hpwl(cells, nets[n]);
                    if (after + 1.0e-6 < before) {
                        std::swap(ids[i], ids[i + 1]);
                        ++pass_accepted;
                    } else {
                        cells[a].x = old_ax;
                        cells[b].x = old_bx;
                    }
                }
            }
        }
        accepted += pass_accepted;
        std::printf("[Legalizer:DREAMPlace] detailed pass=%d accepted_swaps=%zu\n",
                    pass + 1, pass_accepted);
        if (pass_accepted == 0) break;
    }
    return accepted;
}

}  // namespace

LegalizationResult dreamplace_legalize_placement(
    std::vector<Cell>& cells, const std::vector<Net>& nets,
    float_t /*chip_xl*/, float_t /*chip_yl*/, float_t /*chip_xh*/,
    float_t /*chip_yh*/, const std::string& scl_file,
    const LegalizerConfig& config) {
    LegalizationResult result;
    result.method = "dreamplace-greedy";
    for (const Cell& cell : cells) {
        if (!cell.is_terminal) ++result.movable_cells;
    }
    const std::vector<Cell> initial = cells;
    std::vector<Row> rows;
    if (!parse_rows(scl_file, rows, result.message)) return result;
    result.row_spans = rows.size();
    const auto segments = build_blanks(cells, rows);
    for (const auto& row_segments : segments) {
        result.free_segments += row_segments.size();
    }
    if (result.free_segments == 0) {
        result.message = "fixed obstacles leave no legal row segments";
        return result;
    }

    std::printf("[Legalizer:DREAMPlace] Greedy RTL->LTR, alpha=%.3f\n",
                config.greedy_alpha);
    if (!dreamplace_greedy(cells, rows, segments, config.greedy_alpha,
                           result.message) ||
        !check_legality(cells, rows, segments, result.message)) {
        cells = initial;
        return result;
    }
    result.greedy_hpwl = nets.empty() ? 0.0 : compute_hpwl_only(cells, nets);
    const std::vector<Cell> greedy = cells;

    if (config.mode == LegalizerMode::DreamplaceGreedyAbacus) {
        std::string abacus_error;
        if (run_abacus(cells, initial, rows, segments, abacus_error) &&
            check_legality(cells, rows, segments, abacus_error)) {
            result.abacus_hpwl = nets.empty() ? 0.0 : compute_hpwl_only(cells, nets);
            // Unlike a rigid transplant, retain the legal Greedy point when
            // displacement minimization increases the actual nonsmooth HPWL.
            if (!nets.empty() && result.abacus_hpwl > result.greedy_hpwl) {
                cells = greedy;
                std::printf("[Legalizer:DREAMPlace] Abacus HPWL %.1f > Greedy %.1f; "
                            "retaining Greedy\n", result.abacus_hpwl,
                            result.greedy_hpwl);
            } else {
                result.used_abacus = true;
                result.method = "dreamplace-greedy+abacus";
            }
        } else {
            cells = greedy;
            std::printf("[Legalizer:DREAMPlace] Abacus rejected (%s); retaining Greedy\n",
                        abacus_error.c_str());
        }
    }

    result.detailed_swaps = detailed_adjacent_reorder(
        cells, nets, rows, segments, config.detailed_passes);
    if (result.detailed_swaps > 0) {
        std::string detailed_error;
        if (!check_legality(cells, rows, segments, detailed_error)) {
            cells = greedy;
            result.detailed_swaps = 0;
            std::printf("[Legalizer:DREAMPlace] detailed reorder rejected (%s)\n",
                        detailed_error.c_str());
        } else {
            result.detailed_hpwl = compute_hpwl_only(cells, nets);
            result.method += "+exact-reorder";
        }
    }

    displacement(cells, initial, result);
    result.success = true;
    result.message = "ok";
    const double mean = result.movable_cells
        ? result.total_displacement / result.movable_cells : 0.0;
    std::printf("[Legalizer:DREAMPlace] method=%s greedy_hpwl=%.1f "
                "abacus_hpwl=%.1f detailed_hpwl=%.1f swaps=%zu "
                "mean_disp=%.3f max_disp=%.3f\n",
                result.method.c_str(), result.greedy_hpwl, result.abacus_hpwl,
                result.detailed_hpwl, result.detailed_swaps, mean,
                result.max_displacement);
    return result;
}
