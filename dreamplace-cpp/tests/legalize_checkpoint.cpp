#include "bookshelf.h"
#include "legalizer.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace dpcpp;

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) {
        std::cerr << "usage: legalize_checkpoint <raw-base> <generated-global.pl> "
                     "<output.pl> [--dreamplace-detailed]\n";
        return 1;
    }
    if (argc == 5 && std::string(argv[4]) != "--dreamplace-detailed") {
        throw std::runtime_error("unknown legalization option");
    }
    std::string lower = argv[2];
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower.find("eplace") != std::string::npos || lower.find("nsp") != std::string::npos) {
        throw std::runtime_error("forbidden external placement checkpoint");
    }
    Database db = read_bookshelf(argv[1]);
    std::ifstream input(argv[2]);
    if (!input) throw std::runtime_error("cannot open generated checkpoint");
    std::string line;
    int loaded = 0;
    while (std::getline(input, line)) {
        std::istringstream stream(line);
        std::string name;
        Real left = 0.0, bottom = 0.0;
        if (!(stream >> name >> left >> bottom)) continue;
        const auto found = db.node_by_name.find(name);
        if (found == db.node_by_name.end()) continue;
        Node& node = db.nodes[found->second];
        if (!node.fixed) {
            node.x = left + 0.5 * node.width;
            node.y = bottom + 0.5 * node.height;
            ++loaded;
        }
    }
    if (loaded != static_cast<int>(db.movable_ids.size())) {
        throw std::runtime_error("checkpoint does not contain every movable node");
    }
    LegalizeConfig config;
    config.run_dreamplace_detailed = argc == 5;
    const LegalizeResult result = legalize_and_refine(db, config);
    write_bookshelf_pl(db, argv[3]);
    std::cout << std::setprecision(12)
              << "checkpoint_legal_hpwl=" << result.hpwl_after_detailed << '\n';
    return result.legality.legal ? 0 : 2;
}
