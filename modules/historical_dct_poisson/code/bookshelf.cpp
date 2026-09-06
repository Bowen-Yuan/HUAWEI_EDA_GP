#include "bookshelf.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

namespace dpcpp {
namespace {

std::vector<std::string> split(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> fields;
    std::string field;
    while (stream >> field) fields.push_back(field);
    return fields;
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

void reject_forbidden_placement(const std::filesystem::path& path) {
    const std::string name = lower(path.filename().string());
    if (name.find("eplace") != std::string::npos ||
        name.find("nsp") != std::string::npos ||
        name.find("gp.pl") != std::string::npos) {
        throw std::runtime_error("forbidden warm-start placement path: " + path.string());
    }
}

Real value_after_colon(const std::string& line) {
    const auto colon = line.find(':');
    if (colon == std::string::npos) throw std::runtime_error("missing ':' in SCL line");
    return std::stod(line.substr(colon + 1));
}

}  // namespace

std::string basename_of(const std::string& path) {
    return std::filesystem::path(path).filename().string();
}

Database read_bookshelf(const std::string& benchmark_base) {
    namespace fs = std::filesystem;
    Database db;
    db.benchmark_base = fs::absolute(benchmark_base).lexically_normal().string();
    const fs::path base(db.benchmark_base);
    const fs::path aux_path = base.string() + ".aux";
    const fs::path nodes_path = base.string() + ".nodes";
    const fs::path nets_path = base.string() + ".nets";
    const fs::path pl_path = base.string() + ".pl";
    const fs::path scl_path = base.string() + ".scl";
    reject_forbidden_placement(pl_path);
    for (const fs::path& required : {aux_path, nodes_path, nets_path, pl_path, scl_path}) {
        if (!fs::is_regular_file(required)) {
            throw std::runtime_error("missing raw Bookshelf file: " + required.string());
        }
    }
    db.raw_pl_path = fs::absolute(pl_path).lexically_normal().string();

    std::ifstream aux(aux_path);
    std::string line;
    std::getline(aux, line);
    const std::string aux_lower = lower(line);
    if (aux_lower.find("eplace") != std::string::npos ||
        aux_lower.find("nsp") != std::string::npos ||
        aux_lower.find("gp.pl") != std::string::npos) {
        throw std::runtime_error("AUX references a forbidden placement artifact");
    }

    {
        std::ifstream input(nodes_path);
        while (std::getline(input, line)) {
            const auto fields = split(line);
            if (fields.size() < 3 || fields[0] == "UCLA" ||
                fields[0] == "NumNodes" || fields[0] == "NumTerminals" ||
                (!fields[0].empty() && fields[0][0] == '#')) continue;
            Node node;
            node.id = static_cast<int>(db.nodes.size());
            node.name = fields[0];
            node.width = std::stod(fields[1]);
            node.height = std::stod(fields[2]);
            if (fields.size() >= 4) {
                node.fixed = fields[3].rfind("terminal", 0) == 0;
                node.terminal_ni = fields[3] == "terminal_NI";
            }
            db.node_by_name.emplace(node.name, node.id);
            db.nodes.push_back(std::move(node));
        }
    }
    if (db.nodes.empty()) throw std::runtime_error("no nodes parsed");

    {
        std::ifstream input(pl_path);
        std::vector<unsigned char> seen(db.nodes.size(), 0);
        while (std::getline(input, line)) {
            const auto fields = split(line);
            if (fields.size() < 3 || fields[0] == "UCLA" ||
                (!fields[0].empty() && fields[0][0] == '#')) continue;
            const auto found = db.node_by_name.find(fields[0]);
            if (found == db.node_by_name.end()) continue;
            Node& node = db.nodes[found->second];
            node.x = std::stod(fields[1]) + 0.5 * node.width;
            node.y = std::stod(fields[2]) + 0.5 * node.height;
            const auto colon = std::find(fields.begin(), fields.end(), ":");
            if (colon != fields.end() && colon + 1 != fields.end()) node.orientation = *(colon + 1);
            seen[node.id] = 1;
        }
        for (const Node& node : db.nodes) {
            if (node.fixed && !seen[node.id]) {
                throw std::runtime_error("fixed node missing from raw PL: " + node.name);
            }
        }
    }

    {
        std::ifstream input(scl_path);
        Row current;
        bool in_row = false;
        while (std::getline(input, line)) {
            if (line.find("CoreRow") != std::string::npos) {
                current = Row{};
                in_row = true;
            } else if (in_row && line.find("Coordinate") != std::string::npos) {
                current.y = value_after_colon(line);
            } else if (in_row && line.find("Height") != std::string::npos) {
                current.height = value_after_colon(line);
            } else if (in_row && line.find("Sitewidth") != std::string::npos) {
                current.site_width = value_after_colon(line);
            } else if (in_row && line.find("Sitespacing") != std::string::npos) {
                current.site_spacing = value_after_colon(line);
            } else if (in_row && line.find("SubrowOrigin") != std::string::npos) {
                const auto fields = split(line);
                if (fields.size() < 6) throw std::runtime_error("malformed SubrowOrigin");
                Row subrow = current;
                subrow.origin = std::stod(fields[2]);
                subrow.num_sites = std::stoi(fields[5]);
                db.rows.push_back(subrow);
            } else if (in_row && line.find("End") != std::string::npos) {
                in_row = false;
            }
        }
    }
    if (db.rows.empty()) throw std::runtime_error("no rows parsed");
    std::sort(db.rows.begin(), db.rows.end(), [](const Row& a, const Row& b) {
        return a.y < b.y || (a.y == b.y && a.origin < b.origin);
    });
    db.xl = std::numeric_limits<Real>::infinity();
    db.yl = std::numeric_limits<Real>::infinity();
    db.xh = -std::numeric_limits<Real>::infinity();
    db.yh = -std::numeric_limits<Real>::infinity();
    for (const Row& row : db.rows) {
        db.xl = std::min(db.xl, row.origin);
        db.yl = std::min(db.yl, row.y);
        db.xh = std::max(db.xh, row.xh());
        db.yh = std::max(db.yh, row.y + row.height);
    }

    {
        std::ifstream input(nets_path);
        db.node_pin_weight.assign(db.nodes.size(), 0);
        while (std::getline(input, line)) {
            if (line.find("NetDegree") == std::string::npos) continue;
            const auto header = split(line);
            if (header.size() < 3) throw std::runtime_error("malformed NetDegree");
            const int degree = std::stoi(header[2]);
            Net net;
            net.id = static_cast<int>(db.nets.size());
            if (header.size() >= 4) net.name = header[3];
            for (int p = 0; p < degree; ++p) {
                if (!std::getline(input, line)) throw std::runtime_error("truncated net");
                const auto fields = split(line);
                if (fields.empty()) continue;
                const auto found = db.node_by_name.find(fields[0]);
                if (found == db.node_by_name.end()) continue;
                Pin pin;
                pin.node = found->second;
                const auto colon = std::find(fields.begin(), fields.end(), ":");
                if (colon != fields.end() && std::distance(colon, fields.end()) >= 3) {
                    pin.offset_x = std::stod(*(colon + 1));
                    pin.offset_y = std::stod(*(colon + 2));
                }
                net.pins.push_back(pin);
                ++db.node_pin_weight[pin.node];
            }
            if (net.pins.size() >= 2) db.nets.push_back(std::move(net));
        }
    }

    for (const Node& node : db.nodes) {
        (node.fixed ? db.fixed_ids : db.movable_ids).push_back(node.id);
    }
    std::cout << "[Bookshelf] nodes=" << db.nodes.size()
              << " movable=" << db.movable_ids.size()
              << " fixed=" << db.fixed_ids.size()
              << " nets=" << db.nets.size()
              << " rows=" << db.rows.size() << '\n';
    std::cout << "[Bookshelf] raw_pl=" << db.raw_pl_path << '\n';
    std::cout << "[Bookshelf] bounds=[" << db.xl << ',' << db.xh << "]x["
              << db.yl << ',' << db.yh << "]\n";
    return db;
}

void load_movable_placement(Database& db, const std::string& path) {
    namespace fs = std::filesystem;
    const fs::path pl_path = fs::absolute(path).lexically_normal();
    reject_forbidden_placement(pl_path);
    if (!fs::is_regular_file(pl_path)) {
        throw std::runtime_error("missing warm-start placement: " + pl_path.string());
    }

    std::ifstream input(pl_path);
    if (!input) throw std::runtime_error("cannot read " + pl_path.string());
    std::vector<unsigned char> seen(db.nodes.size(), 0);
    std::string line;
    while (std::getline(input, line)) {
        const auto fields = split(line);
        if (fields.size() < 3 || fields[0] == "UCLA" ||
            (!fields[0].empty() && fields[0][0] == '#')) continue;
        const auto found = db.node_by_name.find(fields[0]);
        if (found == db.node_by_name.end()) continue;
        Node& node = db.nodes[found->second];
        if (node.fixed) continue;
        const Real x = std::stod(fields[1]) + 0.5 * node.width;
        const Real y = std::stod(fields[2]) + 0.5 * node.height;
        if (!std::isfinite(x) || !std::isfinite(y) ||
            x < db.xl + 0.5 * node.width || x > db.xh - 0.5 * node.width ||
            y < db.yl + 0.5 * node.height || y > db.yh - 0.5 * node.height) {
            throw std::runtime_error("warm-start node outside placement region: " + node.name);
        }
        node.x = x;
        node.y = y;
        const auto colon = std::find(fields.begin(), fields.end(), ":");
        if (colon != fields.end() && colon + 1 != fields.end()) {
            node.orientation = *(colon + 1);
        }
        seen[node.id] = 1;
    }
    for (int id : db.movable_ids) {
        if (!seen[id]) {
            throw std::runtime_error("movable node missing from warm-start PL: " +
                                     db.nodes[id].name);
        }
    }
    std::cout << "[Init] movable warm start=" << pl_path.string()
              << " nodes=" << db.movable_ids.size() << '\n';
}

void center_gaussian_initialize(Database& db, std::uint64_t seed,
                                Real sigma_ratio) {
    std::mt19937_64 generator(seed);
    const Real cx = 0.5 * (db.xl + db.xh);
    const Real cy = 0.5 * (db.yl + db.yh);
    std::normal_distribution<Real> nx(cx, sigma_ratio * (db.xh - db.xl));
    std::normal_distribution<Real> ny(cy, sigma_ratio * (db.yh - db.yl));
    Real mean_x = 0.0;
    Real mean_y = 0.0;
    for (int id : db.movable_ids) {
        Node& node = db.nodes[id];
        node.x = std::clamp(nx(generator), db.xl + 0.5 * node.width,
                            db.xh - 0.5 * node.width);
        node.y = std::clamp(ny(generator), db.yl + 0.5 * node.height,
                            db.yh - 0.5 * node.height);
        mean_x += node.x;
        mean_y += node.y;
    }
    mean_x /= std::max<std::size_t>(1, db.movable_ids.size());
    mean_y /= std::max<std::size_t>(1, db.movable_ids.size());
    std::cout << "[Init] center Gaussian seed=" << seed
              << " sigma_ratio=" << sigma_ratio
              << " centroid=(" << mean_x << ',' << mean_y << ")\n";
}

void write_bookshelf_pl(const Database& db, const std::string& path) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot write " + path);
    output << "UCLA pl 1.0\n\n";
    output << std::fixed << std::setprecision(6);
    for (const Node& node : db.nodes) {
        output << node.name << ' ' << node.x - 0.5 * node.width << ' '
               << node.y - 0.5 * node.height << " : " << node.orientation;
        if (node.terminal_ni) output << " /FIXED_NI";
        else if (node.fixed) output << " /FIXED";
        output << '\n';
    }
}

}  // namespace dpcpp
