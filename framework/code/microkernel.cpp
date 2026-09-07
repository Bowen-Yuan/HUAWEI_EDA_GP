#include "microkernel.hpp"

#include "epsilon_active/density.hpp"
#include "epsilon_active/hpwl.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#ifdef _WIN32
#include <windows.h>
#include <wincrypt.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;
namespace nsgp {

void ModuleRegistry::add(std::string name, StageFunction function) {
    if (name.empty() || !function || modules_.count(name)) {
        throw std::invalid_argument("invalid or duplicate module registration: " + name);
    }
    modules_.emplace(std::move(name), std::move(function));
}

const StageFunction& ModuleRegistry::get(const std::string& name) const {
    const auto found = modules_.find(name);
    if (found == modules_.end()) throw std::invalid_argument("unknown module: " + name);
    return found->second;
}

std::vector<std::string> ModuleRegistry::names() const {
    std::vector<std::string> result;
    for (const auto& item : modules_) result.push_back(item.first);
    return result;
}

ExactMetrics exact_audit(ea::Database& db, const DensityConfig& config) {
    ea::ExactHpwl hpwl(db);
    ea::ExactOverlapDensity density(
        db, config.bins_x, config.bins_y, config.target_density);
    const auto d = density.evaluate(0.0, 1.0, nullptr, nullptr);
    return {hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr),
            d.energy, d.overflow, d.max_density};
}

void clamp_movable(ea::Database& db) {
    for (const int id : db.movable_ids) {
        auto& node = db.nodes[id];
        node.x = std::clamp(node.x, db.xl + .5 * node.width, db.xh - .5 * node.width);
        node.y = std::clamp(node.y, db.yl + .5 * node.height, db.yh - .5 * node.height);
    }
}

void configure_threads(int threads) {
    if (threads < 1 || threads > 40) throw std::invalid_argument("threads must be 1..40");
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#endif
}

std::string sha256_file(const fs::path& path) {
#ifdef _WIN32
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot hash file: " + path.string());
    HCRYPTPROV provider=0; HCRYPTHASH hash=0;
    if (!CryptAcquireContext(&provider,nullptr,nullptr,PROV_RSA_AES,CRYPT_VERIFYCONTEXT) ||
        !CryptCreateHash(provider,CALG_SHA_256,0,0,&hash)) {
        if (provider) CryptReleaseContext(provider,0);
        throw std::runtime_error("SHA-256 initialization failed");
    }
    std::array<char,1<<15> buffer{};
    while (in.read(buffer.data(),buffer.size()) || in.gcount()) {
        if (!CryptHashData(hash,reinterpret_cast<const BYTE*>(buffer.data()),
                           static_cast<DWORD>(in.gcount()),0)) {
            CryptDestroyHash(hash); CryptReleaseContext(provider,0);
            throw std::runtime_error("SHA-256 update failed");
        }
    }
    std::array<BYTE,32> digest{}; DWORD size=static_cast<DWORD>(digest.size());
    if (!CryptGetHashParam(hash,HP_HASHVAL,digest.data(),&size,0)) {
        CryptDestroyHash(hash); CryptReleaseContext(provider,0);
        throw std::runtime_error("SHA-256 finalization failed");
    }
    CryptDestroyHash(hash); CryptReleaseContext(provider,0);
    std::ostringstream out; out<<std::hex<<std::setfill('0');
    for (DWORD i=0;i<size;++i) out<<std::setw(2)<<static_cast<unsigned>(digest[i]);
    return out.str();
#else
    throw std::runtime_error("SHA-256 is not implemented on this platform");
#endif
}

ExperimentLog::ExperimentLog(const fs::path& output_root,
                             const std::string& run_id,
                             const Json& parameters)
    : root_(fs::absolute(output_root / run_id)) {
    if (run_id.empty()) throw std::invalid_argument("run_id must not be empty");
    if (fs::exists(root_)) throw std::runtime_error("experiment already exists: " + root_.string());
    fs::create_directories(root_);
    std::ofstream(root_ / "params.json") << parameters.dump(2) << '\n';
    trajectory_.open(root_ / "trajectory.csv");
    if (!trajectory_) throw std::runtime_error("cannot create trajectory.csv");
    trajectory_ << "stage_index,module,hpwl_before,hpwl_after,"
                   "overflow_percent_before,overflow_percent_after,"
                   "density_energy,max_density,iterations,accepted,rejected,"
                   "objective_evaluations,wall_seconds\n";
}

void ExperimentLog::record(const StageRecord& r) {
    records_.push_back(r);
    trajectory_ << r.stage_index << ',' << r.module << ',' << std::setprecision(14)
                << r.before.hpwl << ',' << r.after.hpwl << ','
                << r.before.overflow_ratio * 100.0 << ','
                << r.after.overflow_ratio * 100.0 << ','
                << r.after.density_energy << ',' << r.after.max_density << ','
                << r.stats.iterations << ',' << r.stats.accepted << ','
                << r.stats.rejected << ',' << r.stats.objective_evaluations << ','
                << r.wall_seconds << '\n';
    trajectory_.flush();
}

void ExperimentLog::finish(const ExactMetrics& initial,
                           const ExactMetrics& final, double wall_seconds,
                           bool placement_retained) {
    std::ofstream out(root_ / "experiment.md");
    if (!out) throw std::runtime_error("cannot create experiment.md");
    int accepted = 0, rejected = 0, evaluations = 0;
    out << "# Pipeline experiment\n\n## Module chain\n\n";
    for (std::size_t i = 0; i < records_.size(); ++i) {
        if (i) out << " → ";
        out << records_[i].module;
        accepted += records_[i].stats.accepted;
        rejected += records_[i].stats.rejected;
        evaluations += records_[i].stats.objective_evaluations;
    }
    out << "\n\n## Exact metrics\n\n"
        << std::setprecision(14)
        << "- initial HPWL: " << initial.hpwl << "\n"
        << "- initial overflow: " << initial.overflow_ratio * 100.0 << "%\n"
        << "- final HPWL: " << final.hpwl << "\n"
        << "- final overflow: " << final.overflow_ratio * 100.0 << "%\n"
        << "- final density energy: " << final.density_energy << "\n"
        << "- final max density: " << final.max_density << "\n"
        << "\n## Runtime and search\n\n"
        << "- wall seconds: " << wall_seconds << "\n"
        << "- accepted: " << accepted << "\n"
        << "- rejected: " << rejected << "\n"
        << "- objective evaluations: " << evaluations << "\n"
        << "\n## Retention\n\n"
        << (placement_retained
            ? "A placement was retained only by explicit --save-stage request.\n"
            : "Metrics-only: no placement or snapshot retained.\n");
}

void ExperimentLog::fail(const ExactMetrics& initial,
                         const ExactMetrics& last_audited,
                         double wall_seconds, const std::string& message) {
    trajectory_.flush();
    std::ofstream out(root_ / "experiment.md");
    if (!out) throw std::runtime_error("cannot create failure experiment.md");
    out << "# Failed pipeline experiment\n\n"
        << "- status: failed\n"
        << "- error: `" << message << "`\n"
        << std::setprecision(14)
        << "- initial HPWL: " << initial.hpwl << "\n"
        << "- initial overflow: " << initial.overflow_ratio * 100.0 << "%\n"
        << "- last audited HPWL: " << last_audited.hpwl << "\n"
        << "- last audited overflow: " << last_audited.overflow_ratio * 100.0 << "%\n"
        << "- completed stage count: " << records_.size() << "\n"
        << "- wall seconds: " << wall_seconds << "\n\n"
        << "Metrics-only failure record; no placement or snapshot retained.\n";
}

ModuleRegistry make_default_registry() {
    ModuleRegistry registry;
    modules::register_layout_init(registry);
    modules::register_hpwl_adam(registry);
    modules::register_exact_joint_gp(registry);
    modules::register_exact_recovery(registry);
    modules::register_surplus_bisection(registry);
    modules::register_equal_shape_swap(registry);
    modules::register_global_capacity_transport(registry);
    modules::register_density_coordinate(registry);
    modules::register_global_view_gp(registry);
    return registry;
}

}  // namespace nsgp
