#include "epsilon_active/bookshelf.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

namespace ea {
namespace {

std::vector<std::string> fields_of(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> fields;
    for (std::string field; stream >> field;) fields.push_back(field);
    return fields;
}

Real after_colon(const std::string& line) {
    const std::size_t colon = line.find(':');
    if (colon == std::string::npos) {
        throw std::runtime_error("malformed Bookshelf SCL field: " + line);
    }
    return std::stod(line.substr(colon + 1));
}

std::filesystem::path normalize_base(std::filesystem::path path) {
    if (path.extension() == ".aux") path.replace_extension();
    return std::filesystem::absolute(path).lexically_normal();
}

void require_file(const std::filesystem::path& path) {
    if (!std::filesystem::is_regular_file(path)) {
        throw std::runtime_error("missing Bookshelf file: " + path.string());
    }
}

}  // namespace

void load_bookshelf_placement(Database& db, const std::filesystem::path& path) {
    require_file(path);
    std::ifstream input(path);
    std::vector<unsigned char> seen(db.nodes.size(), 0);
    for (std::string line; std::getline(input, line);) {
        const auto fields = fields_of(line);
        if (fields.size() < 3 || fields[0] == "UCLA" || fields[0][0] == '#') {
            continue;
        }
        const auto found = db.node_by_name.find(fields[0]);
        if (found == db.node_by_name.end()) continue;
        Node& node = db.nodes[found->second];
        node.x = std::stod(fields[1]) + 0.5 * node.width;
        node.y = std::stod(fields[2]) + 0.5 * node.height;
        const auto colon = std::find(fields.begin(), fields.end(), ":");
        if (colon != fields.end() && colon + 1 != fields.end()) {
            node.orientation = *(colon + 1);
        }
        seen[node.id] = 1;
    }
    for (const Node& node : db.nodes) {
        if (node.fixed && !seen[node.id]) {
            throw std::runtime_error("fixed node missing from PL: " + node.name);
        }
    }
}

Database read_bookshelf(const std::filesystem::path& benchmark) {
    Database db;
    const std::filesystem::path base = normalize_base(benchmark);
    const std::filesystem::path aux_path = base.string() + ".aux";
    const std::filesystem::path nodes_path = base.string() + ".nodes";
    const std::filesystem::path nets_path = base.string() + ".nets";
    const std::filesystem::path pl_path = base.string() + ".pl";
    const std::filesystem::path scl_path = base.string() + ".scl";
    for (const auto& path : {aux_path, nodes_path, nets_path, pl_path, scl_path}) {
        require_file(path);
    }
    db.benchmark_base = base.string();
    db.raw_pl_path = pl_path.string();

    std::string line;
    {
        std::ifstream input(nodes_path);
        while (std::getline(input, line)) {
            const auto fields = fields_of(line);
            if (fields.size() < 3 || fields[0] == "UCLA" ||
                fields[0] == "NumNodes" || fields[0] == "NumTerminals" ||
                fields[0][0] == '#') {
                continue;
            }
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
    if (db.nodes.empty()) throw std::runtime_error("Bookshelf nodes file is empty");

    load_bookshelf_placement(db, pl_path);

    {
        std::ifstream input(scl_path);
        Row row;
        bool in_row = false;
        while (std::getline(input, line)) {
            if (line.find("CoreRow") != std::string::npos) {
                row = Row{};
                in_row = true;
            } else if (in_row && line.find("Coordinate") != std::string::npos) {
                row.y = after_colon(line);
            } else if (in_row && line.find("Height") != std::string::npos) {
                row.height = after_colon(line);
            } else if (in_row && line.find("Sitespacing") != std::string::npos) {
                row.site_spacing = after_colon(line);
            } else if (in_row && line.find("SubrowOrigin") != std::string::npos) {
                const auto fields = fields_of(line);
                if (fields.size() < 6) {
                    throw std::runtime_error("malformed SubrowOrigin line");
                }
                Row subrow = row;
                subrow.origin = std::stod(fields[2]);
                subrow.num_sites = std::stoi(fields[5]);
                db.rows.push_back(subrow);
            } else if (in_row && line.find("End") != std::string::npos) {
                in_row = false;
            }
        }
    }
    if (db.rows.empty()) throw std::runtime_error("Bookshelf SCL contains no rows");
    db.xl = std::numeric_limits<Real>::infinity();
    db.yl = std::numeric_limits<Real>::infinity();
    db.xh = -std::numeric_limits<Real>::infinity();
    db.yh = -std::numeric_limits<Real>::infinity();
    for (const Row& row : db.rows) {
        db.xl = std::min(db.xl, row.origin);
        db.yl = std::min(db.yl, row.y);
        db.xh = std::max(db.xh, row.right());
        db.yh = std::max(db.yh, row.y + row.height);
    }

    db.node_pin_count.assign(db.nodes.size(), 0);
    {
        std::ifstream input(nets_path);
        while (std::getline(input, line)) {
            if (line.find("NetDegree") == std::string::npos) continue;
            const auto header = fields_of(line);
            if (header.size() < 3) throw std::runtime_error("malformed NetDegree line");
            const int declared_degree = std::stoi(header[2]);
            Net net;
            net.id = static_cast<int>(db.nets.size());
            net.name = header.size() >= 4 ? header[3] : "n" + std::to_string(net.id);
            net.pin_begin = db.pins.size();
            for (int p = 0; p < declared_degree; ++p) {
                if (!std::getline(input, line)) {
                    throw std::runtime_error("truncated Bookshelf net");
                }
                const auto fields = fields_of(line);
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
                db.pins.push_back(pin);
            }
            net.pin_count = db.pins.size() - net.pin_begin;
            if (net.pin_count >= 2) {
                db.nets.push_back(std::move(net));
            } else {
                db.pins.resize(net.pin_begin);
            }
        }
    }

    const std::filesystem::path weights_path = base.string() + ".wts";
    if (std::filesystem::is_regular_file(weights_path)) {
        std::unordered_map<std::string, Real> weights;
        std::ifstream input(weights_path);
        while (std::getline(input, line)) {
            const auto fields = fields_of(line);
            if (fields.size() < 2 || fields[0] == "UCLA" || fields[0][0] == '#') continue;
            try {
                weights[fields[0]] = std::stod(fields[1]);
            } catch (const std::exception&) {
            }
        }
        for (Net& net : db.nets) {
            const auto found = weights.find(net.name);
            if (found != weights.end()) net.weight = found->second;
        }
    }

    for (const Node& node : db.nodes) {
        if (node.fixed) {
            db.fixed_ids.push_back(node.id);
        } else {
            db.movable_ids.push_back(node.id);
            db.movable_area += node.area();
        }
    }

    std::fill(db.node_pin_count.begin(), db.node_pin_count.end(), 0);
    for (const Pin& pin : db.pins) ++db.node_pin_count[pin.node];
    db.node_pin_offsets.resize(db.nodes.size() + 1, 0);
    for (std::size_t i = 0; i < db.nodes.size(); ++i) {
        db.node_pin_offsets[i + 1] = db.node_pin_offsets[i] + db.node_pin_count[i];
    }
    db.node_pin_indices.resize(db.pins.size());
    std::vector<std::size_t> cursor = db.node_pin_offsets;
    for (std::size_t pin = 0; pin < db.pins.size(); ++pin) {
        db.node_pin_indices[cursor[db.pins[pin].node]++] = pin;
    }
    return db;
}

void initialize_center_gaussian(Database& db, std::uint64_t seed,
                                Real sigma_ratio) {
    if (sigma_ratio < 0.0) throw std::invalid_argument("sigma ratio must be nonnegative");
    std::mt19937_64 generator(seed);
    std::normal_distribution<Real> normal(0.0, 1.0);
    const Real center_x = 0.5 * (db.xl + db.xh);
    const Real center_y = 0.5 * (db.yl + db.yh);
    const Real sigma_x = sigma_ratio * (db.xh - db.xl);
    const Real sigma_y = sigma_ratio * (db.yh - db.yl);
    for (int id : db.movable_ids) {
        Node& node = db.nodes[id];
        node.x = std::clamp(center_x + sigma_x * normal(generator),
                            db.xl + 0.5 * node.width,
                            db.xh - 0.5 * node.width);
        node.y = std::clamp(center_y + sigma_y * normal(generator),
                            db.yl + 0.5 * node.height,
                            db.yh - 0.5 * node.height);
    }
}

void write_bookshelf_placement(const Database& db,
                               const std::filesystem::path& path) {
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot write placement: " + path.string());
    output << "UCLA pl 1.0\n# epsilon-active exact nonsmooth global placement\n\n";
    output << std::setprecision(12);
    for (const Node& node : db.nodes) {
        output << node.name << '\t' << node.x - 0.5 * node.width << '\t'
               << node.y - 0.5 * node.height << "\t: " << node.orientation;
        if (node.terminal_ni) output << " /FIXED_NI";
        else if (node.fixed) output << " /FIXED";
        output << '\n';
    }
}

}  // namespace ea
