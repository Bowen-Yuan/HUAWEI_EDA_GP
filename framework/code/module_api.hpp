#pragma once
#include "layout.hpp"
#include <filesystem>
#include <optional>
#include <string>
namespace nsgp {
struct StateArtifact { std::string kind; std::filesystem::path path; std::string file_sha256; };
struct RunContext { std::string run_id, case_name; int stage_index=0, threads=1; std::uint64_t seed=0; std::filesystem::path stage_result_dir; bool save_trajectory=false; int snapshot_every=0; };
struct ModuleStats { int iterations=0, accepted_moves=0, rejected_moves=0; double wall_seconds=0; };
struct ModuleResult { Layout selected_layout,last_layout; bool selected_is_last=true; ModuleStats stats; std::optional<StateArtifact> state_out; std::string status="ok",message; };
}
