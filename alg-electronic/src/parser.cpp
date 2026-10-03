#include "types.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <unordered_map>
#include <cassert>
#include <algorithm>

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
        if (tokens.size() >= 6) {
            nodes_file = bench_name.substr(0, bench_name.find_last_of('/') + 1) + tokens[2];
            nets_file  = bench_name.substr(0, bench_name.find_last_of('/') + 1) + tokens[3];
            pl_file    = bench_name.substr(0, bench_name.find_last_of('/') + 1) + tokens[5];
            scl_file   = bench_name.substr(0, bench_name.find_last_of('/') + 1) + tokens[6];
        } else if (tokens.size() >= 5) {
            // Alternative format
            std::string dir = bench_name.substr(0, bench_name.find_last_of('/') + 1);
            nodes_file = dir + tokens[2];
            nets_file  = dir + tokens[3];
            pl_file    = dir + tokens[5];
            scl_file   = dir + tokens[6];
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
                    // Bookshelf terminals are identified by the explicit
                    // keyword; in ISPD 2005 they are at the end, not first.
                    c.is_terminal = (tok.size() >= 4 &&
                                     tok[3].rfind("terminal", 0) == 0);
                    cells.push_back(c);
                }
            }
        }
        num_terminals_out = (int)std::count_if(
            cells.begin(), cells.end(), [](const Cell& c) { return c.is_terminal; });
        if (num_terminals_out != num_terminals) {
            std::cerr << "Warning: NumTerminals declares " << num_terminals
                      << " but " << num_terminals_out
                      << " terminal records were parsed" << std::endl;
        }
        std::cout << "Parsed " << cells.size() << " nodes ("
                  << num_terminals << " terminals, "
                  << (cells.size() - num_terminals) << " movable)" << std::endl;
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
        std::vector<bool> placed(cells.size(), false);
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
                    // Bookshelf .pl coordinates are lower-left corners;
                    // internal placement coordinates are cell centers.
                    cells[id].x = std::stof(tok[1]) + 0.5f * cells[id].width;
                    cells[id].y = std::stof(tok[2]) + 0.5f * cells[id].height;
                    placed[id] = true;
                }
            }
        }

        // For movable cells still at (0,0), initialize randomly within chip area
        float_t cw = chip_xh - chip_xl;
        float_t ch = chip_yh - chip_yl;
        int init_count = 0;
        for (auto& c : cells) {
            if (!c.is_terminal && !placed[c.id]) {
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
