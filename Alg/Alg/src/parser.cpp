#include "types.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <unordered_map>
#include <cassert>

// Fast string split for tab/space-delimited files
static inline void fast_split(const std::string& line, std::vector<std::string>& out) {
    out.clear();
    const char* s = line.c_str();
    while (*s) {
        while (*s == ' ' || *s == '\t' || *s == '\r') ++s;
        if (!*s) break;
        const char* start = s;
        while (*s && *s != ' ' && *s != '\t' && *s != '\r') ++s;
        out.emplace_back(start, s - start);
    }
}

bool load_placement_file(const std::string& placement_file,
                         std::vector<Cell>& cells)
{
    std::ifstream fin(placement_file);
    if (!fin) {
        std::cerr << "Cannot open placement file: " << placement_file << std::endl;
        return false;
    }

    std::unordered_map<std::string, int> name_to_id;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        name_to_id["o" + std::to_string(i)] = i;
    }
    std::vector<bool> assigned(cells.size(), false);
    std::string line;
    int assigned_count = 0;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#' || line.find("UCLA") != std::string::npos) {
            continue;
        }
        std::vector<std::string> tok;
        fast_split(line, tok);
        if (tok.size() < 3) continue;
        auto it = name_to_id.find(tok[0]);
        if (it == name_to_id.end()) {
            std::cerr << "Unknown node in placement file: " << tok[0] << std::endl;
            return false;
        }
        int id = it->second;
        // Bookshelf .pl uses lower-left coordinates; all Alg solvers and the
        // legalizer use cell centres internally.
        cells[id].x = std::stof(tok[1]) + cells[id].width * 0.5f;
        cells[id].y = std::stof(tok[2]) + cells[id].height * 0.5f;
        if (!assigned[id]) {
            assigned[id] = true;
            ++assigned_count;
        }
    }
    if (assigned_count != static_cast<int>(cells.size())) {
        std::cerr << "Incomplete placement file: assigned " << assigned_count
                  << " of " << cells.size() << " nodes" << std::endl;
        return false;
    }
    std::cout << "Loaded external initial placement: " << placement_file << std::endl;
    return true;
}

bool parse_ispd2005(const std::string& bench_name,
                     std::vector<Cell>& cells,
                     std::vector<Net>& nets,
                     float_t& chip_xl, float_t& chip_yl,
                     float_t& chip_xh, float_t& chip_yh,
                     int& num_terminals_out)
{
    std::string aux_file = bench_name + ".aux";
    std::ifstream aux(aux_file);
    if (!aux.is_open()) {
        std::cerr << "Error: Cannot open " << aux_file << std::endl;
        return false;
    }

    // Parse .aux to get file names
    std::string line, nodes_file, nets_file, wts_file, pl_file, scl_file;
    std::getline(aux, line);
    {
        std::vector<std::string> tokens;
        fast_split(line, tokens);
        // Format: RowBasedPlacement : nodes nets wts pl scl
        if (tokens.size() >= 7) {
            nodes_file = bench_name.substr(0, bench_name.find_last_of("/\\") + 1) + tokens[2];
            nets_file  = bench_name.substr(0, bench_name.find_last_of("/\\") + 1) + tokens[3];
            pl_file    = bench_name.substr(0, bench_name.find_last_of("/\\") + 1) + tokens[5];
            scl_file   = bench_name.substr(0, bench_name.find_last_of("/\\") + 1) + tokens[6];
        } else {
            std::cerr << "Malformed RowBasedPlacement declaration in " << aux_file << std::endl;
            return false;
        }
    }

    // ---- Parse .nodes ----
    {
        std::ifstream fin(nodes_file);
        if (!fin) { std::cerr << "Cannot open " << nodes_file << std::endl; return false; }
        int num_nodes = 0, num_terminals = 0;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (line.find("NumNodes") != std::string::npos) {
                sscanf(line.c_str(), "NumNodes : %d", &num_nodes);
            } else if (line.find("NumTerminals") != std::string::npos) {
                sscanf(line.c_str(), "NumTerminals : %d", &num_terminals);
            } else if (line.find("UCLA") != std::string::npos) {
                continue;
            } else {
                // Data line: name width height [terminal]
                std::vector<std::string> tok;
                fast_split(line, tok);
                if (tok.size() >= 3) {
                    Cell c;
                    c.id = (int)cells.size();
                    c.width = std::stof(tok[1]);
                    c.height = std::stof(tok[2]);
                    c.area = c.width * c.height;
                    c.x = 0; c.y = 0;
                    // Bookshelf does not require terminals to appear first.  The
                    // per-node marker is authoritative; NumTerminals is only a
                    // header count used for validation below.
                    c.is_terminal = tok.size() >= 4 && tok[3].find("terminal") == 0;
                    cells.push_back(c);
                }
            }
        }
        int parsed_terminals = 0;
        for (const auto& c : cells) parsed_terminals += c.is_terminal ? 1 : 0;
        if (parsed_terminals != num_terminals) {
            std::cerr << "Warning: NumTerminals=" << num_terminals
                      << " but parsed " << parsed_terminals << " terminal markers" << std::endl;
        }
        num_terminals_out = parsed_terminals;
        std::cout << "Parsed " << cells.size() << " nodes ("
                  << parsed_terminals << " terminals, "
                  << (cells.size() - parsed_terminals) << " movable)" << std::endl;
    }

    // ---- Parse .scl ----
    {
        std::ifstream fin(scl_file);
        if (!fin) { std::cerr << "Cannot open " << scl_file << std::endl; return false; }

        struct Row { float_t coord, height, sitewidth, sitespacing, origin; int numsites; };
        std::vector<Row> rows;
        Row cur_row;
        bool in_row = false;

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
                    cur_row.coord = std::stof(line.substr(line.find(':') + 1));
                else if (line.find("Height") != std::string::npos)
                    cur_row.height = std::stof(line.substr(line.find(':') + 1));
                else if (line.find("Sitewidth") != std::string::npos)
                    cur_row.sitewidth = std::stof(line.substr(line.find(':') + 1));
                else if (line.find("Sitespacing") != std::string::npos)
                    cur_row.sitespacing = std::stof(line.substr(line.find(':') + 1));
                else if (line.find("SubrowOrigin") != std::string::npos) {
                    // Format: SubrowOrigin : <origin> NumSites : <numsites>
                    // Note: two colons in the line, so tokens: [SubrowOrigin, :, origin, NumSites, :, numsites]
                    std::vector<std::string> tok;
                    fast_split(line, tok);
                    cur_row.origin = std::stof(tok[2]);
                    cur_row.numsites = std::stoi(tok[5]);
                }
            }
        }

        if (rows.empty()) {
            std::cerr << "No rows found in scl file" << std::endl;
            return false;
        }

        chip_yl = rows.front().coord;
        chip_yh = rows.back().coord + rows.back().height;
        chip_xl = rows[0].origin;
        float_t max_x = 0;
        for (auto& r : rows) {
            float_t rx = r.origin + r.numsites * r.sitewidth;
            if (rx > max_x) max_x = rx;
            chip_xl = std::min(chip_xl, r.origin);
        }
        chip_xh = max_x;

        std::cout << "Chip area: [" << chip_xl << ", " << chip_xh << "] x ["
                  << chip_yl << ", " << chip_yh << "], rows=" << rows.size() << std::endl;
    }

    // ---- Parse .pl (try eplace-ip.pl first, then fall back to .pl) ----
    {
        std::string ip_pl_file = bench_name + ".eplace-ip.pl";
        std::ifstream fin(ip_pl_file);
        if (!fin) {
            // Fall back to default .pl
            fin.open(pl_file);
            std::cout << "Using default .pl: " << pl_file << std::endl;
        } else {
            std::cout << "Using ePlace initial placement: " << ip_pl_file << std::endl;
        }

        if (!fin) { std::cerr << "Cannot open placement file" << std::endl; return false; }
        std::unordered_map<std::string, int> name_to_id;
        // Build name mapping
        for (int i = 0; i < (int)cells.size(); ++i) {
            name_to_id["o" + std::to_string(i)] = i;
        }

        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (line.find("UCLA") != std::string::npos) continue;

            std::vector<std::string> tok;
            fast_split(line, tok);
            if (tok.size() >= 3) {
                auto it = name_to_id.find(tok[0]);
                if (it != name_to_id.end()) {
                    int id = it->second;
                    // Convert Bookshelf lower-left to Alg's centre convention.
                    cells[id].x = std::stof(tok[1]) + cells[id].width * 0.5f;
                    cells[id].y = std::stof(tok[2]) + cells[id].height * 0.5f;
                }
            }
        }

        // For movable cells still at (0,0), initialize randomly within chip area
        float_t cw = chip_xh - chip_xl;
        float_t ch = chip_yh - chip_yl;
        int init_count = 0;
        for (auto& c : cells) {
            if (!c.is_terminal && c.x == 0 && c.y == 0) {
                c.x = chip_xl + (float_t(rand()) / RAND_MAX) * cw;
                c.y = chip_yl + (float_t(rand()) / RAND_MAX) * ch;
                ++init_count;
            }
        }
        std::cout << "Initialized " << init_count << " cells with random positions" << std::endl;
    }

    // ---- Parse .nets ----
    {
        std::ifstream fin(nets_file);
        if (!fin) { std::cerr << "Cannot open " << nets_file << std::endl; return false; }

        // Build name-to-id map
        std::unordered_map<std::string, int> name_to_id;
        for (int i = 0; i < (int)cells.size(); ++i) {
            name_to_id["o" + std::to_string(i)] = i;
        }

        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (line.find("UCLA") != std::string::npos) continue;
            if (line.find("NumNets") != std::string::npos) continue;
            if (line.find("NumPins") != std::string::npos) continue;

            if (line.find("NetDegree") != std::string::npos) {
                Net net;
                net.id = (int)nets.size();

                // Parse NetDegree:  K  netname
                std::vector<std::string> tok;
                fast_split(line, tok);
                int degree = std::stoi(tok[2]);

                // Read degree pin lines
                for (int p = 0; p < degree; ++p) {
                    if (!std::getline(fin, line)) break;
                    // Format: cell_name  I/O : x_offset y_offset
                    std::vector<std::string> ptok;
                    fast_split(line, ptok);
                    if (ptok.size() >= 1) {
                        auto it = name_to_id.find(ptok[0]);
                        if (it != name_to_id.end()) {
                            net.cell_ids.push_back(it->second);
                            // <cell> <I/O> : <x-offset> <y-offset>
                            net.pin_offset_x.push_back(ptok.size() >= 5 ? std::stof(ptok[3]) : 0.0f);
                            net.pin_offset_y.push_back(ptok.size() >= 5 ? std::stof(ptok[4]) : 0.0f);
                        }
                    }
                }
                if (!net.cell_ids.empty()) {
                    nets.push_back(std::move(net));
                }
            }
        }
        std::cout << "Parsed " << nets.size() << " nets" << std::endl;
    }

    return true;
}
