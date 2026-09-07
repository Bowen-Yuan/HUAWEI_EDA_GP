#include "global_view_lab.hpp"
#include "microkernel.hpp"

#include "epsilon_active/bookshelf.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/hpwl.hpp"
#include "epsilon_active/optimizer.hpp"
#include <json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Json = nlohmann::json;
#ifndef NSGP_GIT_COMMIT
#define NSGP_GIT_COMMIT "unknown"
#endif
#ifndef NSGP_GIT_BRANCH
#define NSGP_GIT_BRANCH "unknown"
#endif
namespace {
constexpr const char* kDefaultDataset = "D:\\codex_project\\HUAWEI_EDA\\alg-electronic\\ispd2005";
constexpr const char* kDefaultCheckpoint =
    "D:\\codex_project\\HUAWEI_EDA\\epsilon-active\\experiments\\h375_a1_surplus_recovery10_exchange\\run\\global.pl";

struct LabOptions {
    std::string case_name = "adaptec1";
    fs::path dataset = kDefaultDataset;
    fs::path checkpoint = kDefaultCheckpoint;
    fs::path output_root = "framework/results/experiments";
    std::string run_id;
    std::string optimizer = "adam";
    std::string step_policy = "control";
    int iterations = 50, threads = 1, bins_x = 512, bins_y = 512;
    double target_density = 1.0, learning_rate = 0.02, lambda = 0.20;
    double max_delta_bins = 0.50, trust_radius_bins = 1.0, bundle_mix = .5;
    bool bundle = false, active_ensemble = false, capacity_transport = false;
    double transport_target_percent = 6.0;
    std::optional<std::string> save_stage;
};

struct Metrics { double hpwl=0, energy=0, overflow=0, max_density=0; };

std::string now_id() {
    const auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm tm{}; localtime_s(&tm, &t); char text[32];
    std::strftime(text, sizeof(text), "%Y%m%d_%H%M%S", &tm); return text;
}

LabOptions parse(int argc, char** argv) {
    LabOptions o;
    for (int i=2; i<argc; ++i) {
        const std::string a=argv[i]; auto value=[&]() { if (++i>=argc) throw std::runtime_error("missing value for "+a); return std::string(argv[i]); };
        if (a=="--case") o.case_name=value(); else if (a=="--dataset-root") o.dataset=value();
        else if (a=="--input-checkpoint") o.checkpoint=value(); else if (a=="--output-root") o.output_root=value();
        else if (a=="--run-id") o.run_id=value(); else if (a=="--optimizer") o.optimizer=value();
        else if (a=="--step-policy") o.step_policy=value(); else if (a=="--iterations") o.iterations=std::stoi(value());
        else if (a=="--threads") o.threads=std::stoi(value()); else if (a=="--lambda") o.lambda=std::stod(value());
        else if (a=="--learning-rate") o.learning_rate=std::stod(value()); else if (a=="--bundle") o.bundle=true;
        else if (a=="--active-ensemble") o.active_ensemble=true; else if (a=="--bundle-mix") o.bundle_mix=std::stod(value());
        else if (a=="--capacity-transport") o.capacity_transport=true;
        else if (a=="--transport-target-percent") o.transport_target_percent=std::stod(value());
        else if (a=="--save-stage") o.save_stage=value(); else throw std::runtime_error("unknown lab option: "+a);
    }
    if (o.run_id.empty()) o.run_id=now_id()+"_"+o.case_name+"_"+o.optimizer;
    if (o.iterations<1 || o.threads<1 || o.bins_x<1 || o.bins_y<1 || o.target_density<=0 || o.target_density>1)
        throw std::runtime_error("invalid global-view lab parameters");
    if (o.step_policy!="control" && o.step_policy!="trust") throw std::runtime_error("step policy must be control or trust");
    return o;
}

class ExperimentWorkspace {
public:
    explicit ExperimentWorkspace(std::string id) {
        root_ = fs::temp_directory_path()/"nonsmooth-gp"/id;
        fs::create_directories(root_);
        const auto normalized = root_.lexically_normal();
        if (normalized.parent_path().filename() != "nonsmooth-gp") throw std::runtime_error("unsafe experiment workspace");
    }
    ~ExperimentWorkspace() { std::error_code ec; fs::remove_all(root_, ec); }
    const fs::path& root() const { return root_; }
private: fs::path root_;
};

Metrics evaluate(ea::Database& db, int bx, int by, double target) {
    ea::ExactHpwl hpwl(db); ea::ExactOverlapDensity density(db,bx,by,target);
    const auto dm=density.evaluate(0.0,1.0,nullptr,nullptr);
    return {hpwl.evaluate(0.0,1.0,-1,nullptr,nullptr),dm.energy,dm.overflow,dm.max_density};
}
double dot(const std::vector<double>& a,const std::vector<double>& b) { double v=0; for(size_t i=0;i<a.size();++i)v+=a[i]*b[i]; return v; }
void normalize(std::vector<double>& x) { const double n=std::sqrt(std::max(1e-30,dot(x,x)/std::max<size_t>(1,x.size()))); for(double&v:x)v/=n; }
std::vector<double> simplex_project(std::vector<double> v) {
    auto u=v; std::sort(u.rbegin(),u.rend()); double total=0, theta=0; int rho=0;
    for(size_t j=0;j<u.size();++j) { total+=u[j]; const double q=(total-1.0)/(j+1); if(u[j]>q){rho=int(j)+1;theta=q;} }
    if(!rho) return std::vector<double>(v.size(),1.0/v.size());
    for(double& x:v) x=std::max(0.0,x-theta);
    return v;
}
std::vector<double> min_norm(const std::vector<std::vector<double>>& dirs) {
    const size_t k=dirs.size(); if(k==1) return dirs[0]; std::vector<double> a(k,1.0/k), q(k*k);
    double lipschitz=1e-12; for(size_t i=0;i<k;++i) for(size_t j=0;j<k;++j) { q[i*k+j]=dot(dirs[i],dirs[j]); lipschitz=std::max(lipschitz,std::abs(q[i*k+j])); }
    for(int step=0;step<80;++step) { std::vector<double> g(k); for(size_t i=0;i<k;++i)for(size_t j=0;j<k;++j)g[i]+=q[i*k+j]*a[j]; for(size_t i=0;i<k;++i)a[i]-=g[i]/(lipschitz*k); a=simplex_project(a); }
    std::vector<double> out(dirs[0].size()); for(size_t i=0;i<k;++i)for(size_t j=0;j<out.size();++j)out[j]+=a[i]*dirs[i][j]; return out;
}
void clamp(ea::Database& db) { for(int id:db.movable_ids) { auto& n=db.nodes[id]; n.x=std::clamp(n.x,db.xl+.5*n.width,db.xh-.5*n.width); n.y=std::clamp(n.y,db.yl+.5*n.height,db.yh-.5*n.height); } }
void write_placement(const ea::Database& db,const fs::path& out) { std::ofstream f(out); f<<"UCLA pl 1.0\n\n"<<std::setprecision(12); for(const auto& n:db.nodes) f<<n.name<<'\t'<<n.x-.5*n.width<<'\t'<<n.y-.5*n.height<<"\t: "<<n.orientation<<(n.fixed?" /FIXED":"")<<'\n'; }

void write_metadata(const fs::path& root,const LabOptions& o,const std::string& input_hash,const Metrics& initial) {
    Json chain=Json::array();
    if(o.capacity_transport) chain.push_back({{"module","global_capacity_transport"},{"target_percent",o.transport_target_percent},{"candidate_grid",64},{"audit_grid",o.bins_x}});
    chain.push_back({{"module","global_view_gp"},{"active_ensemble",o.active_ensemble},{"epsilon_bin_scales",o.active_ensemble?Json{0,.25,.5,1.0}:Json{0}},{"temporal_bundle",o.bundle},{"bundle_size",o.bundle?4:0},{"bundle_mix",o.bundle_mix}});
    Json p={{"run_id",o.run_id},{"case",o.case_name},{"threads",o.threads},{"seed",219},{"git",{{"branch",NSGP_GIT_BRANCH},{"commit",NSGP_GIT_COMMIT},{"dirty","not_checked"}}},{"input_checkpoint",fs::absolute(o.checkpoint).string()},{"input_checkpoint_sha256",input_hash},{"initial_exact_metrics",{{"hpwl",initial.hpwl},{"overflow_percent",initial.overflow*100},{"density_energy",initial.energy},{"max_density",initial.max_density}}},{"density_grid",{{"bins_x",o.bins_x},{"bins_y",o.bins_y},{"target_density",o.target_density}}},{"module_chain",chain},{"optimizer",{{"name",o.optimizer},{"learning_rate",o.learning_rate},{"maximum_delta_bins",o.max_delta_bins}}},{"step_policy",{{"name",o.step_policy},{"trust_radius_bins",o.trust_radius_bins},{"max_backtracks",9}}},{"acceptance",{{"policy","strict_exact_cap_and_hpwl_decrease"},{"overflow_cap_percent",initial.overflow*100+1e-5},{"lambda",o.lambda}}},{"iterations",o.iterations},{"retention",{{"keep_final_placement",false},{"keep_snapshots",false},{"keep_debug_artifacts",false},{"explicit_saved_stage",o.save_stage.value_or("")},{"cleanup_workspace_on_success",true},{"cleanup_workspace_on_failure",true}}}};
    std::ofstream(root/"params.json")<<p.dump(2)<<'\n';
}

void write_summary(const fs::path& root,const LabOptions& o,const std::string& hash,
                   const Metrics& initial,const Metrics& last,const Metrics& best,
                   int transport_moves,int accepted,int rejected,int resets,double seconds) {
    std::ofstream md(root/"experiment.md");
    md<<"# V3 global-view experiment\n\n## Provenance\n\n"
      <<"- input checkpoint: `"<<fs::absolute(o.checkpoint).string()<<"`\n"
      <<"- input SHA-256: `"<<hash<<"`\n- threads: "<<o.threads<<"\n"
      <<"\n## Module chain\n\n";
    if(o.capacity_transport) md<<"global_capacity_transport → ";
    md<<"global_view_gp\n\n## Exact metrics\n\n"<<std::setprecision(14)
      <<"- initial: HPWL "<<initial.hpwl<<", overflow "<<initial.overflow*100<<"%\n"
      <<"- best feasible: HPWL "<<best.hpwl<<", overflow "<<best.overflow*100<<"%\n"
      <<"- last: HPWL "<<last.hpwl<<", overflow "<<last.overflow*100<<"%\n"
      <<"\n## Runtime and search\n\n- wall seconds: "<<seconds
      <<"\n- capacity moves: "<<transport_moves<<"\n- accepted steps: "<<accepted
      <<"\n- rejected steps: "<<rejected<<"\n- optimizer/bundle resets: "<<resets
      <<"\n- NaN/exception: no\n\n## Retention\n\n"
      <<(o.save_stage?"One stage placement explicitly retained.\n":"Metrics-only; no placement or snapshot retained.\n");
}

// Coarse grid is used only to nominate a legal move.  Every candidate is
// audited on the canonical grid before acceptance; macros are excluded.
int capacity_transport(ea::Database& db, const LabOptions& o, Metrics& current) {
    constexpr int kCoarse = 64; const double target=o.transport_target_percent/100.0;
    ea::ExactOverlapDensity coarse(db,kCoarse,kCoarse,o.target_density);
    const double bw=(db.xh-db.xl)/kCoarse, bh=(db.yh-db.yl)/kCoarse; int accepted=0;
    for(int pass=0; pass<24 && current.overflow>target; ++pass) {
        coarse.evaluate(0,1,nullptr,nullptr); const auto& occ=coarse.occupancy();
        int source=-1,dest=-1; double excess=0,best=std::numeric_limits<double>::infinity();
        const double cap=coarse.bin_area()*o.target_density;
        for(int i=0;i<kCoarse*kCoarse;++i) { if(occ[i]-cap>excess){excess=occ[i]-cap;source=i;} if(occ[i]<best){best=occ[i];dest=i;} }
        if(source<0 || dest<0 || source==dest) break;
        const int sx=source%kCoarse,sy=source/kCoarse,dx=dest%kCoarse,dy=dest/kCoarse;
        int chosen=-1; for(int id:db.movable_ids) { const auto& n=db.nodes[id]; if(n.width>bw || n.height>bh) continue; const int x=std::clamp(int((n.x-db.xl)/bw),0,kCoarse-1), y=std::clamp(int((n.y-db.yl)/bh),0,kCoarse-1); if(x==sx&&y==sy){chosen=id;break;} }
        if(chosen<0) break;
        auto& n=db.nodes[chosen];
        const double ox=n.x,oy=n.y;
        n.x=db.xl+(dx+.5)*bw; n.y=db.yl+(dy+.5)*bh; clamp(db); Metrics candidate=evaluate(db,o.bins_x,o.bins_y,o.target_density);
        if(candidate.overflow<current.overflow && candidate.hpwl<=current.hpwl*1.002) { current=candidate; ++accepted; }
        else { n.x=ox;n.y=oy; }
    }
    return accepted;
}
}

int run_global_view_lab(int argc, char** argv) {
    try {
        const LabOptions o=parse(argc,argv); nsgp::configure_threads(o.threads); if(!fs::is_regular_file(o.checkpoint)) throw std::runtime_error("external input checkpoint does not exist: "+o.checkpoint.string());
        const std::string input_hash=nsgp::sha256_file(o.checkpoint); ExperimentWorkspace workspace(o.run_id);
        ea::Database db=ea::read_bookshelf(o.dataset/o.case_name/o.case_name); ea::load_bookshelf_placement(db,o.checkpoint); clamp(db);
        const Metrics initial=evaluate(db,o.bins_x,o.bins_y,o.target_density); const double cap=initial.overflow+1e-7;
        const fs::path root=fs::absolute(o.output_root/o.run_id); if(fs::exists(root)) throw std::runtime_error("experiment result directory already exists: "+root.string()); fs::create_directories(root);
        write_metadata(root,o,input_hash,initial); std::ofstream tr(root/"trajectory.csv");
        tr<<"iteration,hpwl,overflow_percent,max_density,density_energy,best_feasible_hpwl,best_feasible_overflow_percent,lambda,step,trust_radius,accepted,wall_seconds,direction_cosine,bundle_ratio\n";
        ea::ExactHpwl hpwl(db); ea::ExactOverlapDensity density(db,o.bins_x,o.bins_y,o.target_density);
        auto opt=ea::make_optimizer(ea::parse_optimizer(o.optimizer),.9,.999,.9,1e-8); const size_t n=db.nodes.size(); opt->reset(2*n);
        Metrics current=initial,best=initial; const int transport_moves=o.capacity_transport?capacity_transport(db,o,current):0; best=current; double radius=o.trust_radius_bins*std::min(density.bin_width(),density.bin_height()); int reject_streak=0, accept_streak=0,total_accepted=0,total_rejected=0,total_resets=0; std::vector<std::vector<double>> history;
        const auto start=std::chrono::steady_clock::now();
        for(int it=1;it<=o.iterations;++it) {
            std::vector<std::vector<double>> wires; const double bin=std::min(density.bin_width(),density.bin_height());
            const std::vector<double> scales=o.active_ensemble?std::vector<double>{0,.25,.5,1.0}:std::vector<double>{0};
            for(double s:scales) { std::vector<double> gx,gy; hpwl.evaluate(s*bin,1.0,-1,&gx,&gy); std::vector<double> d(2*n); for(size_t j=0;j<n;++j){d[j]=gx[j];d[n+j]=gy[j];} normalize(d); wires.push_back(std::move(d)); }
            std::vector<double> wire=min_norm(wires), dx,dy; density.evaluate(0.0,1.0,&dx,&dy); std::vector<double> g(2*n); for(size_t j=0;j<n;++j){g[j]=wire[j]+o.lambda*dx[j];g[n+j]=wire[n+j]+o.lambda*dy[j];} normalize(g);
            double cosine=1.0,bundle_ratio=1.0;
            if(o.bundle) { if(history.size()==4) history.erase(history.begin()); history.push_back(g); auto b=min_norm(history); const double bn=std::sqrt(dot(b,b)),gn=std::sqrt(dot(g,g)); cosine=dot(b,g)/std::max(1e-30,bn*gn); bundle_ratio=bn/std::max(1e-30,gn); for(size_t j=0;j<g.size();++j)g[j]=(1-o.bundle_mix)*g[j]+o.bundle_mix*b[j]; normalize(g); }
            std::vector<double> delta; const double max_delta=o.step_policy=="trust"?radius:o.max_delta_bins*bin; opt->compute_delta(g,o.learning_rate,max_delta,delta);
            std::vector<std::pair<double,double>> old; old.reserve(db.movable_ids.size()); for(int id:db.movable_ids)old.push_back({db.nodes[id].x,db.nodes[id].y});
            Metrics candidate=current; bool accepted=false;
            // Exact backtracking is deliberately evaluated on the canonical
            // grid; it is not a coarse surrogate acceptance filter.
            for(int backtrack=0; backtrack<9 && !accepted; ++backtrack) {
                const double scale=std::ldexp(1.0,-backtrack);
                for(size_t j=0;j<db.movable_ids.size();++j){auto& node=db.nodes[db.movable_ids[j]];node.x=old[j].first-scale*delta[db.movable_ids[j]];node.y=old[j].second-scale*delta[n+db.movable_ids[j]];} clamp(db);
                candidate=evaluate(db,o.bins_x,o.bins_y,o.target_density);
                accepted=std::isfinite(candidate.hpwl)&&std::isfinite(candidate.overflow)&&candidate.overflow<=cap&&candidate.hpwl<current.hpwl;
            }
            if(accepted) { current=candidate; if(candidate.hpwl<best.hpwl)best=candidate; ++accept_streak; ++total_accepted; reject_streak=0; if(o.step_policy=="trust"&&accept_streak>=3){radius*=1.25;accept_streak=0;} }
            else { for(size_t j=0;j<db.movable_ids.size();++j){auto& node=db.nodes[db.movable_ids[j]];node.x=old[j].first;node.y=old[j].second;} ++reject_streak;++total_rejected;accept_streak=0; if(o.step_policy=="trust")radius*=.5; if(reject_streak>=4){history.clear();opt->reset(2*n);reject_streak=0;++total_resets;} }
            const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            tr<<it<<','<<std::setprecision(14)<<current.hpwl<<','<<current.overflow*100<<','<<current.max_density<<','<<current.energy<<','<<best.hpwl<<','<<best.overflow*100<<','<<o.lambda<<','<<o.learning_rate<<','<<radius<<','<<(accepted?1:0)<<','<<wall<<','<<cosine<<','<<bundle_ratio<<'\n';
        }
        if(o.save_stage) { fs::create_directories(root/"saved"); write_placement(db,root/"saved"/(*o.save_stage+".pl")); }
        const double total_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        write_summary(root,o,input_hash,initial,current,best,transport_moves,total_accepted,total_rejected,total_resets,total_seconds);
        std::cout<<std::setprecision(14)<<"completed metrics-only global-view experiment "<<root<<" HPWL="<<best.hpwl<<" overflow_percent="<<best.overflow*100<<"% capacity_moves="<<transport_moves<<" input_sha256="<<input_hash<<'\n'; return 0;
    } catch(const std::exception& e) { std::cerr<<"nsgp lab: "<<e.what()<<'\n'; return 1; }
}
