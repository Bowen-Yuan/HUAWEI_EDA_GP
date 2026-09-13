#include "microkernel.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <vector>
 #include <system_error>

namespace {

std::vector<std::vector<std::string>> read_csv(
    const std::filesystem::path& path) {
    std::ifstream in(path);
    std::vector<std::vector<std::string>> rows;
    for (std::string line; std::getline(in, line);) {
        std::vector<std::string> fields;
        std::string field;
        std::istringstream stream(line);
        while (std::getline(stream, field, ',')) fields.push_back(field);
        // count trailing empty fields dropped by getline
        std::size_t commas = 0;
        for (char c : line) if (c == ',') ++commas;
        while (fields.size() < commas + 1) fields.emplace_back("");
        rows.push_back(fields);
    }
    return rows;
}

}  // namespace

int main() {
    namespace fs = std::filesystem;
    const auto root = fs::temp_directory_path() / "nsgp_explog_contracts";
    std::error_code ec;
    fs::remove_all(root, ec);
    // A stale locked directory from a previous crashed run must not break
    // the contract: fall back to a fresh run id instead.
    const std::string run_id = fs::exists(root / "contract_run")
        ? "contract_run_" + std::to_string(
              std::chrono::steady_clock::now().time_since_epoch().count())
        : "contract_run";
    nsgp::Json params = {{"run_id", "contract"}};
    nsgp::ExperimentLog log(root, run_id, params);
    // two stages x two iterations, with crossing global iterations
    int global = 0;
    for (int stage = 0; stage < 2; ++stage) {
        for (int it = 0; it < 2; ++it) {
            log.record_iteration(stage, "finite_radius_oracle_gp", nsgp::Json{
                {"stage_iteration", it},
                {"global_iteration", global++},
                {"hpwl", 100.0 + it},
                {"overflow_percent", 8.5},
                {"density_energy", 12.0},
                {"max_density", 1.1},
                {"objective_evaluations", 4 + global},
                {"local_density_queries", 64 * global},
                {"local_hpwl_queries", 0}});
        }
        nsgp::StageRecord record;
        record.stage_index = stage;
        record.module = "finite_radius_oracle_gp";
        record.before = {100.0, 12.0, 0.085, 1.1};
        record.after = {99.0, 11.0, 0.075, 1.05};
        record.stats = {2, 0, 0, 8, 128, 0};
        record.wall_seconds = 0.5;
        log.record(record);
    }
    log.finish({100.0, 12.0, 0.085, 1.1}, {99.0, 11.0, 0.075, 1.05}, 1.0);

    const auto rows = read_csv(root / run_id / "trajectory.csv");
    if (rows.size() < 6) return 1;
    const std::size_t columns = rows[0].size();
    // header: 19 columns ending at wall_seconds
    if (columns != 19) return 2;
    if (rows[0].at(15) != "objective_evaluations") return 3;
    if (rows[0].at(16) != "local_density_queries") return 4;
    if (rows[0].at(17) != "local_hpwl_queries") return 5;
    for (std::size_t r = 1; r < rows.size(); ++r) {
        if (rows[r].size() != columns) return 6;
    }
    // iteration row: overflow values must sit in columns 8 and 9
    if (rows[1].at(8).empty() || rows[1].at(9).empty()) return 7;
    if (std::abs(std::stod(rows[1].at(8)) - 8.5) > 1e-12) return 8;
    if (std::abs(std::stod(rows[1].at(9)) - 8.5) > 1e-12) return 9;
    // stage row: hpwl_before at column 6, hpwl_after at column 7
    std::size_t stage_row = 3;
    if (rows[stage_row].at(0) != "stage") return 10;
    if (std::abs(std::stod(rows[stage_row].at(6)) - 100.0) > 1e-9) return 11;
    if (std::abs(std::stod(rows[stage_row].at(7)) - 99.0) > 1e-9) return 12;
    if (std::abs(std::stod(rows[stage_row].at(8)) - 8.5) > 1e-9) return 13;
    if (std::abs(std::stod(rows[stage_row].at(9)) - 7.5) > 1e-9) return 14;
    // global_iteration monotonic across stages
    long long previous = -1;
    for (std::size_t r = 1; r < rows.size(); ++r) {
        if (rows[r].at(0) != "iteration") continue;
        const long long g = std::stoll(rows[r].at(3));
        if (g <= previous) return 15;
        previous = g;
    }
    // local density queries recorded on stage rows
    if (std::stoll(rows[stage_row].at(16)) != 128) return 16;
    fs::remove_all(root, ec);
    std::cout << "experiment log contracts passed\n";
}
