#pragma once

#include "epsilon_active/types.hpp"
#include <json.hpp>

#include <filesystem>
#include <functional>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace nsgp {

using Json = nlohmann::json;

struct DensityConfig {
    int bins_x = 512;
    int bins_y = 512;
    double target_density = 1.0;
};

struct ExactMetrics {
    double hpwl = 0.0;
    double density_energy = 0.0;
    double overflow_ratio = 0.0;
    double max_density = 0.0;
};

struct StageContext {
    ea::Database& db;
    DensityConfig density;
    int threads = 1;
};

struct StageStats {
    int iterations = 0;
    int accepted = 0;
    int rejected = 0;
    int objective_evaluations = 0;
};

using StageFunction = std::function<StageStats(StageContext&, const Json&)>;

class ModuleRegistry {
public:
    void add(std::string name, StageFunction function);
    const StageFunction& get(const std::string& name) const;
    std::vector<std::string> names() const;
private:
    std::map<std::string, StageFunction> modules_;
};

ExactMetrics exact_audit(ea::Database& db, const DensityConfig& density);
void clamp_movable(ea::Database& db);
void configure_threads(int threads);
std::string sha256_file(const std::filesystem::path& path);

struct StageRecord {
    int stage_index = 0;
    std::string module;
    ExactMetrics before;
    ExactMetrics after;
    StageStats stats;
    double wall_seconds = 0.0;
};

class ExperimentLog {
public:
    ExperimentLog(const std::filesystem::path& output_root,
                  const std::string& run_id, const Json& parameters);
    void record(const StageRecord& record);
    void finish(const ExactMetrics& initial, const ExactMetrics& final,
                double wall_seconds, bool placement_retained = false);
    void fail(const ExactMetrics& initial, const ExactMetrics& last_audited,
              double wall_seconds, const std::string& message);
    const std::filesystem::path& root() const noexcept { return root_; }
private:
    std::filesystem::path root_;
    std::ofstream trajectory_;
    std::vector<StageRecord> records_;
};

namespace modules {
void register_layout_init(ModuleRegistry& registry);
void register_hpwl_adam(ModuleRegistry& registry);
void register_exact_joint_gp(ModuleRegistry& registry);
void register_exact_recovery(ModuleRegistry& registry);
void register_surplus_bisection(ModuleRegistry& registry);
void register_equal_shape_swap(ModuleRegistry& registry);
void register_global_capacity_transport(ModuleRegistry& registry);
void register_density_coordinate(ModuleRegistry& registry);
void register_global_view_gp(ModuleRegistry& registry);
}

ModuleRegistry make_default_registry();

}  // namespace nsgp
