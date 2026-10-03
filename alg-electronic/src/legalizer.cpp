#include "legalizer.h"
#include "dreamplace_legalizer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

float_t compute_hpwl_only(const std::vector<Cell>& cells,
                          const std::vector<Net>& nets);

namespace {

constexpr double kEps = 1.0e-5;

struct Subrow {
    double origin = 0.0;
    int num_sites = 0;
};

struct ParsedRow {
    double y = 0.0;
    double height = 0.0;
    double site_width = 0.0;
    double site_spacing = 0.0;
    std::vector<Subrow> subrows;
};

struct RowSpan {
    double y = 0.0;
    double height = 0.0;
    double site_spacing = 1.0;
    double origin = 0.0;
    double xh = 0.0;
};

struct Interval {
    double xl = 0.0;
    double xh = 0.0;
};

struct Segment {
    int row_span = -1;
    double xl = 0.0;
    double xh = 0.0;
    double site_spacing = 1.0;
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

bool parse_rows(const std::string& scl_file, std::vector<RowSpan>& spans,
                std::string& error) {
    std::ifstream fin(scl_file);
    if (!fin) {
        error = "cannot open SCL file: " + scl_file;
        return false;
    }

    ParsedRow row;
    bool in_row = false;
    std::string line;
    while (std::getline(fin, line)) {
        const auto tok = split_ws(line);
        if (tok.empty() || tok[0] == "#" || tok[0] == "UCLA") continue;

        if (tok[0] == "CoreRow") {
            row = ParsedRow{};
            in_row = true;
            continue;
        }
        if (!in_row) continue;
        if (tok[0] == "End") {
            const double spacing = row.site_spacing > 0.0
                ? row.site_spacing : row.site_width;
            if (row.height <= 0.0 || spacing <= 0.0 || row.subrows.empty()) {
                error = "incomplete CoreRow in SCL file";
                return false;
            }
            for (const auto& subrow : row.subrows) {
                if (subrow.num_sites <= 0) continue;
                spans.push_back(RowSpan{
                    row.y,
                    row.height,
                    spacing,
                    subrow.origin,
                    subrow.origin + subrow.num_sites * spacing
                });
            }
            in_row = false;
            continue;
        }

        if (tok[0] == "Coordinate" && tok.size() >= 3) {
            row.y = std::stod(tok[2]);
        } else if (tok[0] == "Height" && tok.size() >= 3) {
            row.height = std::stod(tok[2]);
        } else if (tok[0] == "Sitewidth" && tok.size() >= 3) {
            row.site_width = std::stod(tok[2]);
        } else if (tok[0] == "Sitespacing" && tok.size() >= 3) {
            row.site_spacing = std::stod(tok[2]);
        } else if (tok[0] == "SubrowOrigin" && tok.size() >= 6) {
            row.subrows.push_back(Subrow{std::stod(tok[2]), std::stoi(tok[5])});
        }
    }

    if (spans.empty()) {
        error = "no placement rows found in SCL file";
        return false;
    }
    std::sort(spans.begin(), spans.end(), [](const RowSpan& a, const RowSpan& b) {
        if (a.y != b.y) return a.y < b.y;
        return a.origin < b.origin;
    });
    return true;
}

double snap_up(double x, double origin, double spacing) {
    return origin + std::ceil((x - origin) / spacing - kEps) * spacing;
}

double snap_down(double x, double origin, double spacing) {
    return origin + std::floor((x - origin) / spacing + kEps) * spacing;
}

void build_free_segments(const std::vector<Cell>& cells,
                         const std::vector<RowSpan>& rows,
                         std::vector<Segment>& segments) {
    for (int row_id = 0; row_id < static_cast<int>(rows.size()); ++row_id) {
        const auto& row = rows[row_id];
        std::vector<Interval> blocked;
        for (const auto& cell : cells) {
            if (!cell.is_terminal) continue;
            const double tyl = static_cast<double>(cell.y) - 0.5 * cell.height;
            const double tyh = static_cast<double>(cell.y) + 0.5 * cell.height;
            if (tyl >= row.y + row.height - kEps || tyh <= row.y + kEps) continue;

            const double txl = static_cast<double>(cell.x) - 0.5 * cell.width;
            const double txh = static_cast<double>(cell.x) + 0.5 * cell.width;
            const double xl = std::max(row.origin, txl);
            const double xh = std::min(row.xh, txh);
            if (xh > xl + kEps) blocked.push_back({xl, xh});
        }

        std::sort(blocked.begin(), blocked.end(), [](const Interval& a, const Interval& b) {
            if (a.xl != b.xl) return a.xl < b.xl;
            return a.xh < b.xh;
        });
        std::vector<Interval> merged;
        for (const auto& interval : blocked) {
            if (merged.empty() || interval.xl > merged.back().xh + kEps) {
                merged.push_back(interval);
            } else {
                merged.back().xh = std::max(merged.back().xh, interval.xh);
            }
        }

        double cursor = row.origin;
        auto add_segment = [&](double raw_xl, double raw_xh) {
            const double xl = snap_up(raw_xl, row.origin, row.site_spacing);
            const double xh = snap_down(raw_xh, row.origin, row.site_spacing);
            if (xh > xl + kEps) {
                segments.push_back(Segment{
                    row_id, xl, xh, row.site_spacing
                });
            }
        };
        for (const auto& interval : merged) {
            add_segment(cursor, interval.xl);
            cursor = std::max(cursor, interval.xh);
        }
        add_segment(cursor, row.xh);
    }
}

bool greedy_blank_pass(std::vector<Cell>& cells,
                       const std::vector<Cell>& anchors,
                       const std::vector<RowSpan>& rows,
                       const std::vector<Segment>& initial_segments,
                       bool left_to_right,
                       std::string& error) {
    std::vector<std::vector<Interval>> row_blanks(rows.size());
    for (const auto& segment : initial_segments) {
        row_blanks[segment.row_span].push_back({segment.xl, segment.xh});
    }
    for (auto& blanks : row_blanks) {
        std::sort(blanks.begin(), blanks.end(), [](const Interval& a, const Interval& b) {
            return a.xl < b.xl;
        });
    }

    std::vector<int> movable;
    movable.reserve(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        if (!cells[i].is_terminal) movable.push_back(i);
    }
    std::sort(movable.begin(), movable.end(), [&](int a, int b) {
        if (anchors[a].x != anchors[b].x) {
            return left_to_right ? anchors[a].x < anchors[b].x
                                 : anchors[a].x > anchors[b].x;
        }
        if (anchors[a].width != anchors[b].width) {
            return anchors[a].width > anchors[b].width;
        }
        return a < b;
    });

    for (int cell_id : movable) {
        const auto& anchor = anchors[cell_id];
        const double raw_left = static_cast<double>(anchor.x) - 0.5 * anchor.width;
        const double anchor_y = anchor.y;
        int best_row = -1;
        int best_blank = -1;
        double best_left = 0.0;
        double best_cost = std::numeric_limits<double>::infinity();
        double best_blank_width = -1.0;

        for (int row_id = 0; row_id < static_cast<int>(rows.size()); ++row_id) {
            const auto& row = rows[row_id];
            if (anchor.height > row.height + kEps) continue;
            const double target_y = row.y + 0.5 * anchor.height;
            const double y_cost = std::abs(anchor_y - target_y);
            if (y_cost > best_cost + kEps) continue;

            auto& blanks = row_blanks[row_id];
            for (int blank_id = 0; blank_id < static_cast<int>(blanks.size()); ++blank_id) {
                const auto& blank = blanks[blank_id];
                const double blank_width = blank.xh - blank.xl;
                if (blank_width + kEps < anchor.width) continue;

                double target_left = snap_down(raw_left, blank.xl, row.site_spacing);
                target_left = std::clamp(target_left, blank.xl,
                                         blank.xh - anchor.width);
                const double tolerance = std::min(
                    4.0 * static_cast<double>(anchor.width), blank_width / 4.0);
                if (target_left <= blank.xl + tolerance) {
                    target_left = blank.xl;
                } else if (target_left + anchor.width >= blank.xh - tolerance) {
                    target_left = blank.xh - anchor.width;
                }

                const double cost = y_cost + std::abs(target_left - raw_left);
                if (cost < best_cost - kEps ||
                    (std::abs(cost - best_cost) <= kEps &&
                     blank_width > best_blank_width)) {
                    best_row = row_id;
                    best_blank = blank_id;
                    best_left = target_left;
                    best_cost = cost;
                    best_blank_width = blank_width;
                }
            }
        }

        if (best_row < 0) {
            error = "greedy blank search cannot place cell o" +
                    std::to_string(anchor.id) + " (width=" +
                    std::to_string(anchor.width) + ")";
            return false;
        }

        auto& blanks = row_blanks[best_row];
        const Interval selected = blanks[best_blank];
        const double placed_right = best_left + anchor.width;
        if (best_left <= selected.xl + kEps) {
            if (placed_right >= selected.xh - kEps) {
                blanks.erase(blanks.begin() + best_blank);
            } else {
                blanks[best_blank].xl = placed_right;
            }
        } else if (placed_right >= selected.xh - kEps) {
            blanks[best_blank].xh = best_left;
        } else {
            blanks[best_blank].xh = best_left;
            blanks.insert(blanks.begin() + best_blank + 1,
                          Interval{placed_right, selected.xh});
        }

        cells[cell_id].x = static_cast<float_t>(best_left + 0.5 * anchor.width);
        cells[cell_id].y = static_cast<float_t>(
            rows[best_row].y + 0.5 * anchor.height);
    }
    return true;
}

bool dreamplace_style_greedy(std::vector<Cell>& cells,
                             const std::vector<Net>& nets,
                             const std::vector<RowSpan>& rows,
                             const std::vector<Segment>& initial_segments,
                             std::string& error) {
    const std::vector<Cell> original = cells;

    struct Candidate {
        std::string name;
        std::vector<Cell> placement;
        float_t hpwl = 0.0f;
        double displacement = 0.0;
    };
    std::vector<Candidate> candidates;

    auto add_candidate = [&](const char* name, const std::vector<Cell>& anchor,
                             bool left_to_right) -> bool {
        std::vector<Cell> trial = original;
        std::string pass_error;
        if (!greedy_blank_pass(trial, anchor, rows, initial_segments,
                               left_to_right, pass_error)) {
            if (error.empty()) error = pass_error;
            return false;
        }
        Candidate candidate;
        candidate.name = name;
        candidate.placement.swap(trial);
        candidate.hpwl = nets.empty()
            ? 0.0f : compute_hpwl_only(candidate.placement, nets);
        for (int i = 0; i < static_cast<int>(original.size()); ++i) {
            if (original[i].is_terminal) continue;
            candidate.displacement += std::hypot(
                static_cast<double>(candidate.placement[i].x - original[i].x),
                static_cast<double>(candidate.placement[i].y - original[i].y));
        }
        candidates.push_back(std::move(candidate));
        return true;
    };

    const bool have_right = add_candidate("right-to-left", original, false);
    const bool have_left = add_candidate("left-to-right", original, true);

    const std::size_t base_count = candidates.size();
    for (std::size_t i = 0; i < base_count; ++i) {
        const Candidate& base = candidates[i];
        std::vector<Cell> blended = original;
        for (int j = 0; j < static_cast<int>(blended.size()); ++j) {
            if (blended[j].is_terminal) continue;
            blended[j].x = 0.5f * (original[j].x + base.placement[j].x);
            blended[j].y = 0.5f * (original[j].y + base.placement[j].y);
        }
        const bool opposite_direction = base.name == "right-to-left";
        const std::string name = base.name + "+opposite";
        add_candidate(name.c_str(), blended, opposite_direction);
    }

    if (candidates.empty()) {
        if (error.empty()) error = "all greedy blank sweeps failed";
        return false;
    }
    auto best = std::min_element(candidates.begin(), candidates.end(),
        [&](const Candidate& a, const Candidate& b) {
            if (!nets.empty() && a.hpwl != b.hpwl) return a.hpwl < b.hpwl;
            return a.displacement < b.displacement;
        });
    const std::size_t movable_count = static_cast<std::size_t>(std::count_if(
        original.begin(), original.end(),
        [](const Cell& cell) { return !cell.is_terminal; }));
    for (const auto& candidate : candidates) {
        const double mean = candidate.displacement /
            std::max<std::size_t>(1, movable_count);
        std::printf("[Legalizer] candidate=%s exact_hpwl=%.1f mean_disp=%.3f%s\n",
                    candidate.name.c_str(), candidate.hpwl, mean,
                    &candidate == &*best ? " [selected]" : "");
    }
    cells.swap(best->placement);
    error.clear();
    (void)have_right;
    (void)have_left;
    return true;
}

}  // namespace

LegalizationResult legalize_placement_legacy(std::vector<Cell>& cells,
                                             const std::vector<Net>& nets,
                                             float_t /*chip_xl*/, float_t /*chip_yl*/,
                                             float_t /*chip_xh*/, float_t /*chip_yh*/,
                                             const std::string& scl_file) {
    LegalizationResult result;
    for (const auto& cell : cells) {
        if (!cell.is_terminal) ++result.movable_cells;
    }
    std::printf("[Legalizer] Starting obstacle-aware row legalization...\n");

    std::vector<RowSpan> rows;
    if (!parse_rows(scl_file, rows, result.message)) {
        std::printf("[Legalizer] ERROR: %s\n", result.message.c_str());
        return result;
    }
    result.row_spans = rows.size();

    std::vector<Segment> segments;
    build_free_segments(cells, rows, segments);
    result.free_segments = segments.size();
    if (segments.empty()) {
        result.message = "fixed obstacles leave no legal row segments";
        std::printf("[Legalizer] ERROR: %s\n", result.message.c_str());
        return result;
    }

    const std::vector<Cell> original = cells;
    if (!dreamplace_style_greedy(cells, nets, rows, segments, result.message)) {
        std::printf("[Legalizer] ERROR: %s\n", result.message.c_str());
        return result;
    }
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        if (cells[i].is_terminal) continue;
        const double dx = static_cast<double>(cells[i].x) - original[i].x;
        const double dy = static_cast<double>(cells[i].y) - original[i].y;
        const double displacement = std::hypot(dx, dy);
        result.total_displacement += displacement;
        result.max_displacement = std::max(result.max_displacement, displacement);
    }

    result.success = true;
    result.method = "legacy-blank-sweep";
    result.message = "ok";
    const double mean_displacement = result.movable_cells > 0
        ? result.total_displacement / result.movable_cells : 0.0;
    std::printf("[Legalizer] Legalized %zu cells in %zu free segments; "
                "mean displacement=%.3f, max=%.3f\n",
                result.movable_cells, result.free_segments,
                mean_displacement, result.max_displacement);
    return result;
}

LegalizationResult legalize_placement(std::vector<Cell>& cells,
                                      const std::vector<Net>& nets,
                                      float_t chip_xl, float_t chip_yl,
                                      float_t chip_xh, float_t chip_yh,
                                      const std::string& scl_file) {
    return legalize_placement_legacy(cells, nets, chip_xl, chip_yl,
                                     chip_xh, chip_yh, scl_file);
}

LegalizationResult legalize_placement(std::vector<Cell>& cells,
                                      const std::vector<Net>& nets,
                                      float_t chip_xl, float_t chip_yl,
                                      float_t chip_xh, float_t chip_yh,
                                      const std::string& scl_file,
                                      const LegalizerConfig& config) {
    if (config.mode == LegalizerMode::Legacy) {
        return legalize_placement_legacy(cells, nets, chip_xl, chip_yl,
                                         chip_xh, chip_yh, scl_file);
    }
    const std::vector<Cell> checkpoint = cells;
    LegalizationResult result = dreamplace_legalize_placement(
        cells, nets, chip_xl, chip_yl, chip_xh, chip_yh, scl_file, config);
    if (result.success || !config.fallback_to_legacy) return result;

    const std::string primary_error = result.message;
    cells = checkpoint;
    std::printf("[Legalizer] DREAMPlace path failed (%s); using legacy fallback\n",
                primary_error.c_str());
    result = legalize_placement_legacy(cells, nets, chip_xl, chip_yl,
                                       chip_xh, chip_yh, scl_file);
    result.used_fallback = true;
    if (!result.success) {
        result.message = "DREAMPlace path: " + primary_error +
                         "; legacy path: " + result.message;
    }
    return result;
}
