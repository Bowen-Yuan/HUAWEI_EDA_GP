#pragma once
#include "common.hpp"
#include <string>
#include <vector>
namespace nsgp {
struct NodeInfo { std::string name; Real width=0, height=0; bool movable=false, fixed=false, terminal_ni=false; };
struct Problem { std::string case_name; std::filesystem::path benchmark_base; std::vector<NodeInfo> nodes; std::vector<NodeId> movable_ids, fixed_ids; Rect region; Real movable_area=0; };
}
