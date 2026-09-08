#include "global_view_lab.hpp"
#include "microkernel.hpp"
#include "epsilon_active/bookshelf.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace fs = std::filesystem;
using nsgp::Json;
#ifndef NSGP_GIT_COMMIT
#define NSGP_GIT_COMMIT "unknown"
#endif
#ifndef NSGP_GIT_BRANCH
#define NSGP_GIT_BRANCH "unknown"
#endif
#ifndef NSGP_GIT_DIRTY
#define NSGP_GIT_DIRTY "not_checked"
#endif
#ifndef NSGP_GIT_REMOTE
#define NSGP_GIT_REMOTE "unknown"
#endif
namespace {
constexpr const char* kDefaultDataset =
    "D:\\codex_project\\HUAWEI_EDA\\alg-electronic\\ispd2005";

struct Options {
    std::string command;
    std::string case_name = "adaptec1";
    fs::path pipeline, placement, dataset = kDefaultDataset;
    fs::path output_root = "framework/results/experiments";
    std::string run_id, save_stage;
    int threads = 1;
    bool all_cases = false;
};

void usage() {
    std::cout << "nsgp run|batch|audit|list-modules|lab [--case adaptec1] "
                 "[--pipeline FILE] [--placement FILE] [--dataset-root DIR] "
                 "[--threads N] [--run-id ID] [--output-root DIR] "
                 "[--save-stage MODULE]\n";
}

Options parse_options(int argc, char** argv) {
    if (argc < 2) { usage(); throw std::runtime_error("missing command"); }
    Options o; o.command = argv[1];
    for (int i=2; i<argc; ++i) {
        const std::string arg=argv[i];
        auto value=[&]() { if (++i>=argc) throw std::runtime_error("missing value for "+arg); return std::string(argv[i]); };
        if (arg=="--case") o.case_name=value();
        else if (arg=="--pipeline") o.pipeline=value();
        else if (arg=="--placement" || arg=="--initial-placement") o.placement=value();
        else if (arg=="--dataset-root") o.dataset=value();
        else if (arg=="--threads") o.threads=std::stoi(value());
        else if (arg=="--run-id") o.run_id=value();
        else if (arg=="--output-root") o.output_root=value();
        else if (arg=="--save-stage") o.save_stage=value();
        else if (arg=="--all-cases") o.all_cases=true;
        else throw std::runtime_error("unknown argument: "+arg);
    }
    nsgp::configure_threads(o.threads);
    return o;
}

Json read_json(const fs::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read JSON: "+path.string());
    Json result; in>>result; return result;
}

std::string timestamp() {
    const auto ticks=std::chrono::system_clock::now().time_since_epoch().count();
    return std::to_string(ticks);
}

std::string wall_clock_time() {
    const auto now=std::chrono::system_clock::now();
    const auto value=std::chrono::system_clock::to_time_t(now);
    std::tm local{};
#ifdef _WIN32
    localtime_s(&local,&value);
#else
    localtime_r(&value,&local);
#endif
    std::ostringstream out;
    out<<std::put_time(&local,"%Y-%m-%dT%H:%M:%S%z");
    return out.str();
}

nsgp::DensityConfig density_config(const Json& pipeline) {
    nsgp::DensityConfig result;
    const Json* value=nullptr;
    if (pipeline.contains("density_grid")) value=&pipeline.at("density_grid");
    else if (pipeline.contains("density")) value=&pipeline.at("density");
    if (value) {
        result.bins_x=value->value("bins_x",result.bins_x);
        result.bins_y=value->value("bins_y",result.bins_y);
        result.target_density=value->value("target_density",result.target_density);
    }
    if (result.bins_x<1 || result.bins_y<1 || result.target_density<=0 ||
        result.target_density>1) throw std::runtime_error("invalid density configuration");
    return result;
}

struct PreparedStage { std::string module; fs::path config_path; Json config; };

std::vector<PreparedStage> prepare_stages(const fs::path& pipeline_path,
                                          const Json& pipeline) {
    std::vector<PreparedStage> stages;
    const std::string compliance=pipeline.value("compliance_mode",std::string("challenge_nonsmooth"));
    for (const auto& stage : pipeline.at("stages")) {
        const std::string module=stage.at("module").get<std::string>();
        if (compliance=="challenge_nonsmooth" && module=="historical_dct_poisson")
            throw std::runtime_error("challenge_nonsmooth rejects historical smooth surrogate");
        const fs::path config=fs::weakly_canonical(
            pipeline_path.parent_path()/stage.at("config").get<std::string>());
        stages.push_back({module,config,read_json(config)});
    }
    if (stages.empty()) throw std::runtime_error("pipeline must contain at least one stage");
    return stages;
}

void run_case(const Options& o,const std::string& case_name) {
    if (o.pipeline.empty()) throw std::runtime_error("--pipeline is required");
    const fs::path pipeline_path=fs::absolute(o.pipeline);
    const Json pipeline=read_json(pipeline_path);
    const auto density=density_config(pipeline);
    const auto stages=prepare_stages(pipeline_path,pipeline);
    auto registry=nsgp::make_default_registry();
    if (!o.save_stage.empty() &&
        std::none_of(stages.begin(),stages.end(),[&](const PreparedStage& stage) {
            return stage.module==o.save_stage;
        })) throw std::runtime_error("--save-stage is not present in pipeline: "+o.save_stage);

    ea::Database db=ea::read_bookshelf(o.dataset/case_name/case_name);
    if (!o.placement.empty()) ea::load_bookshelf_placement(db,o.placement);
    nsgp::clamp_movable(db);

    Json recorded_stages=Json::array();
    for (const auto& stage:stages) recorded_stages.push_back({
        {"module",stage.module},{"config_path",stage.config_path.string()},
        {"parameters",stage.config}});
    const std::string run_id=(o.run_id.empty()?case_name+"_"+timestamp():o.run_id)+
        (case_name==o.case_name?"":"_"+case_name);
    Json params={{"run_id",run_id},{"started_at",wall_clock_time()},
        {"case",case_name},
        {"dataset_root",fs::absolute(o.dataset).string()},
        {"pipeline",pipeline_path.string()},{"threads",o.threads},
        {"git",{{"remote",NSGP_GIT_REMOTE},{"branch",NSGP_GIT_BRANCH},{"commit",NSGP_GIT_COMMIT},
                {"dirty",NSGP_GIT_DIRTY}}},
        {"experiment",pipeline.value("experiment",Json::object())},
        {"input",o.placement.empty()?Json{{"kind","raw_bookshelf"}}:
            Json{{"kind","placement"},{"path",fs::absolute(o.placement).string()},
                 {"sha256",nsgp::sha256_file(o.placement)}}},
        {"exact_evaluator",{{"bins_x",density.bins_x},{"bins_y",density.bins_y},
                            {"target_density",density.target_density},
                            {"overflow_unit","percent"}}},
        {"module_chain",recorded_stages},
        {"retention",{{"metrics_only",o.save_stage.empty()},
                      {"explicit_saved_stage",o.save_stage}}}};

    nsgp::ExperimentLog log(o.output_root,run_id,params);
    const auto started=std::chrono::steady_clock::now();
    const auto initial=nsgp::exact_audit(db,density);
    auto last_audited=initial;
    bool retained=false;
    try {
        for (std::size_t index=0; index<stages.size(); ++index) {
            const auto& stage=stages[index];
            const auto before=nsgp::exact_audit(db,density);
            const auto stage_started=std::chrono::steady_clock::now();
            nsgp::StageContext context{db,density,o.threads,
                [&](int stage_iteration, const std::string& module, const Json& telemetry) {
                    Json tagged=telemetry;
                    tagged["stage_iteration"]=stage_iteration;
                    tagged["global_iteration"]=stage_iteration;
                    log.record_iteration(static_cast<int>(index),module,tagged);
                }};
            const auto stats=registry.get(stage.module)(context,stage.config);
            nsgp::clamp_movable(db);
            const auto after=nsgp::exact_audit(db,density);
            last_audited=after;
            const double seconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-stage_started).count();
            log.record({static_cast<int>(index),stage.module,before,after,stats,seconds});
            if (!o.save_stage.empty() && o.save_stage==stage.module) {
                fs::create_directories(log.root()/"saved");
                ea::write_bookshelf_placement(db,log.root()/"saved"/(stage.module+".pl"));
                retained=true;
            }
        }
        const auto final=nsgp::exact_audit(db,density);
        const double seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();
        log.finish(initial,final,seconds,retained);
        std::cout<<std::setprecision(14)<<"completed "<<case_name<<" -> "<<log.root()
                 <<" HPWL="<<final.hpwl<<" overflow_percent="
                 <<final.overflow_ratio*100.0<<"%\n";
    } catch (const std::exception& error) {
        const double seconds=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-started).count();
        log.fail(initial,last_audited,seconds,error.what());
        throw;
    }
}
}

int main(int argc,char** argv) {
    if (argc>=2 && std::string(argv[1])=="lab") return run_global_view_lab(argc,argv);
    try {
        const auto options=parse_options(argc,argv);
        if (options.command=="list-modules") {
            for (const auto& name:nsgp::make_default_registry().names()) std::cout<<name<<'\n';
            std::cout<<"historical_dct_poisson (historical-only)\n";
            return 0;
        }
        if (options.command=="audit") {
            if (options.placement.empty()) throw std::runtime_error("audit requires --placement");
            ea::Database db=ea::read_bookshelf(options.dataset/options.case_name/options.case_name);
            ea::load_bookshelf_placement(db,options.placement);
            const auto m=nsgp::exact_audit(db,{});
            std::cout<<std::setprecision(14)<<"hpwl="<<m.hpwl
                     <<" density_energy="<<m.density_energy
                     <<" overflow_percent="<<m.overflow_ratio*100.0<<"%"
                     <<" max_density="<<m.max_density<<'\n';
            return 0;
        }
        if (options.command=="run") { run_case(options,options.case_name); return 0; }
        if (options.command=="batch" && options.all_cases) {
            for (const auto& name:std::vector<std::string>{"adaptec1","adaptec2","adaptec3","adaptec4","bigblue1","bigblue2","bigblue3","bigblue4"}) run_case(options,name);
            return 0;
        }
        usage(); throw std::runtime_error("batch requires --all-cases");
    } catch (const std::exception& error) {
        std::cerr<<"nsgp: "<<error.what()<<'\n'; return 1;
    }
}
