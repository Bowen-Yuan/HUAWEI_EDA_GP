#include "epsilon_active/bookshelf.hpp"
#include "epsilon_active/bisection.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/hpwl.hpp"
#include "epsilon_active/placer.hpp"
#include "epsilon_active/recovery.hpp"
#include "epsilon_active/swap_recovery.hpp"
#include "experiment_lab.hpp"
#include <json.hpp>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace fs=std::filesystem; using Json=nlohmann::json;
namespace {
constexpr const char* kDefaultRoot="D:\\codex_project\\HUAWEI_EDA\\alg-electronic\\ispd2005";
struct DensityCfg { int bx=64,by=64; double target=.9; };
struct Opt { std::string command; std::string cas="adaptec1"; fs::path pipeline; fs::path placement; fs::path dataset=kDefaultRoot; int threads=1; bool all=false; };
Json read_json(const fs::path& p){ std::ifstream in(p); if(!in) throw std::runtime_error("cannot read JSON: "+p.string()); Json t; in>>t; return t; }
template<class T> T val(const Json& t,const std::string& k,T d){ return t.contains(k) ? t.at(k).get<T>() : d; }
std::string esc(std::string s){ for(char& c:s) if(c=='\\') c='/'; return s; }
void usage(){ std::cout<<"nsgp run|batch|audit|list-modules [--case adaptec1] [--pipeline FILE] [--placement FILE] [--dataset-root DIR] [--threads N]\n"; }
Opt options(int ac,char** av){ if(ac<2){usage(); throw std::runtime_error("missing command");} Opt o; o.command=av[1]; for(int i=2;i<ac;++i){std::string a=av[i]; auto get=[&](){if(++i>=ac)throw std::runtime_error("missing value for "+a);return std::string(av[i]);}; if(a=="--case")o.cas=get(); else if(a=="--pipeline")o.pipeline=get(); else if(a=="--placement"||a=="--initial-placement")o.placement=get(); else if(a=="--dataset-root")o.dataset=get(); else if(a=="--threads")o.threads=std::stoi(get()); else if(a=="--all-cases")o.all=true; else throw std::runtime_error("unknown argument: "+a);} if(o.threads<1||o.threads>40)throw std::runtime_error("threads must be 1..40"); return o; }
DensityCfg density_cfg(const Json& t){ DensityCfg d; const Json* q=nullptr; if(t.contains("density_grid"))q=&t.at("density_grid"); else if(t.contains("density"))q=&t.at("density"); if(q){d.bx=val(*q,"bins_x",d.bx);d.by=val(*q,"bins_y",d.by);d.target=val(*q,"target_density",d.target);} if(d.bx<1||d.by<1||d.target<=0||d.target>1)throw std::runtime_error("invalid density grid"); return d; }
struct M { double hpwl=0,energy=0,overflow=0,maxd=0; };
M audit(ea::Database& db,DensityCfg d){ ea::ExactHpwl h(db); ea::ExactOverlapDensity den(db,d.bx,d.by,d.target); auto dm=den.evaluate(0,1,nullptr,nullptr); return {h.evaluate(0,1,-1,nullptr,nullptr),dm.energy,dm.overflow,dm.max_density}; }
void csv(std::ofstream& f,int i,const std::string& mod,const M& before,const M& after,double secs,int its){f<<i<<','<<mod<<','<<std::setprecision(14)<<before.hpwl<<','<<after.hpwl<<','<<before.overflow*100<<','<<after.overflow*100<<','<<its<<','<<secs<<"\n";}
void write_pl(const ea::Database& db, const fs::path& path) { std::ofstream out(path); if(!out) throw std::runtime_error("cannot write placement: "+path.string()); out<<"UCLA pl 1.0\n# nsgp exact nonsmooth checkpoint\n\n"<<std::setprecision(12); for(const auto& n:db.nodes) { out<<n.name<<'\t'<<n.x-.5*n.width<<'\t'<<n.y-.5*n.height<<"\t: "<<n.orientation; if(n.terminal_ni) out<<" /FIXED_NI"; else if(n.fixed) out<<" /FIXED"; out<<'\n'; } }
void write_stage(const fs::path& dir, ea::Database& db, const M& m, const std::string& mod, int idx, double sec) {
  std::error_code ec; fs::create_directory(dir, ec); if(ec && ec != std::errc::file_exists) throw std::runtime_error("cannot create stage directory: "+ec.message());
  const fs::path selected=dir / "selected.pl"; write_pl(db, selected);
  Json manifest={{"module",mod},{"stage_index",idx},{"runtime_seconds",sec},{"metrics",{{"hpwl",m.hpwl},{"density_energy",m.energy},{"overflow_percent",m.overflow*100},{"max_density",m.maxd}}}};
  std::ofstream(dir / "stage_manifest.json") << manifest.dump(2) << '\n';
}
void clamp(ea::Database& db){for(int id:db.movable_ids){auto& n=db.nodes[id];n.x=std::clamp(n.x,db.xl+.5*n.width,db.xh-.5*n.width);n.y=std::clamp(n.y,db.yl+.5*n.height,db.yh-.5*n.height);}}
int run_module(const std::string& name,const Json& c,ea::Database& db,DensityCfg d,int threads){ auto nested=[&](const char* a,const char* b,double fallback){return c.contains(a)&&c[a].contains(b)?c[a][b].get<double>():fallback;}; if(name=="layout_init"){auto mode=val<std::string>(c,"mode","raw");if(mode=="center_gaussian")ea::initialize_center_gaussian(db,val<std::uint64_t>(c,"seed",219),val<double>(c,"sigma_ratio",.001));else if(mode!="raw")throw std::runtime_error("layout_init mode must be raw or center_gaussian");return 0;} if(name=="hpwl_adam"||name=="exact_joint_gp"){ea::PlaceConfig pc;pc.bins_x=d.bx;pc.bins_y=d.by;pc.target_density=d.target;pc.threads=threads;pc.iterations=val<int>(c,"iterations",10);pc.hpwl_epsilon=nested("hpwl_direction","epsilon",val<double>(c,"epsilon",125));pc.active_power=nested("hpwl_direction","active_power",val<double>(c,"active_power",4));pc.degree_limit=c.contains("hpwl_direction")?val<int>(c["hpwl_direction"],"degree_limit",val<int>(c,"degree_limit",100)):val<int>(c,"degree_limit",100);pc.step_fraction=nested("optimizer","learning_rate",val<double>(c,"learning_rate",.003));pc.max_step_multiplier=nested("optimizer","maximum_delta",val<double>(c,"maximum_delta",4));pc.lambda.density_weight_scale=name=="hpwl_adam"?0.0:nested("lambda","density_weight_scale",.5);pc.stop_overflow=nested("selector","overflow_cap",val<double>(c,"overflow_cap",.07));pc.batch_acceptance.enabled=c.contains("exact_batch_acceptance")?val<bool>(c["exact_batch_acceptance"],"enabled",false):false;pc.batch_acceptance.max_trials=c.contains("exact_batch_acceptance")?val<int>(c["exact_batch_acceptance"],"max_backtracks",0):0;pc.coarse_flow.bins_x=d.bx;pc.coarse_flow.bins_y=d.by;pc.coarse_flow.passes=0;pc.output_dir=fs::path("framework/results")/"_scratch";pc.snapshot_every=0; ea::global_place(db,pc);clamp(db);return pc.iterations;} if(name=="exact_recovery"){ea::ExactOverlapDensity den(db,d.bx,d.by,d.target);ea::RecoveryConfig rc;rc.sweeps=val<int>(c,"sweeps",1);rc.line_search_steps=val<int>(c,"line_search_steps",4);rc.step_bins=val<double>(c,"step_bins",8);rc.overflow_cap=val<double>(c,"overflow_cap",.07);rc.net_block=val<bool>(c,"net_block",false);rc.compact_directions=val<bool>(c,"compact_directions",false);ea::recover_hpwl_under_overflow(db,den,rc);return rc.sweeps;} if(name=="surplus_bisection"){ea::ExactOverlapDensity den(db,d.bx,d.by,d.target);ea::BisectionConfig bc;bc.enabled=true;bc.leaf_bins=val<int>(c,"leaf_bins",8);bc.position_seeded=val<bool>(c,"position_seeded",true);bc.surplus_only=val<bool>(c,"surplus_only",true);bc.leaf_nearest_capacity=val<bool>(c,"leaf_nearest_capacity",true);bc.leaf_hpwl_guided=val<bool>(c,"leaf_hpwl_guided",true);ea::recursive_hypergraph_bisection(db,den,bc);return 1;} if(name=="equal_shape_swap"){ea::ExactOverlapDensity den(db,d.bx,d.by,d.target);ea::SwapRecoveryConfig sc;sc.sweeps=val<int>(c,"sweeps",1);sc.radius_bins=val<int>(c,"radius_bins",4);sc.candidates=val<int>(c,"candidates",16);sc.exact_shortlist=val<int>(c,"exact_shortlist",0);sc.net_aware=val<bool>(c,"net_aware",true);ea::recover_hpwl_with_equal_shape_swaps(db,den,sc);return sc.sweeps;} if(name=="historical_dct_poisson")throw std::runtime_error("historical_dct_poisson is deliberately not built in V1");throw std::runtime_error("unknown module: "+name); }
void run_case(const Opt& o, const std::string& cas) {
  if(o.pipeline.empty()) throw std::runtime_error("--pipeline is required for run/batch");
  auto pipe=read_json(o.pipeline); auto compliance=val<std::string>(pipe,"compliance_mode","challenge_nonsmooth");
  ea::Database db=ea::read_bookshelf(o.dataset/cas/cas); if(!o.placement.empty()) ea::load_bookshelf_placement(db,o.placement);
  DensityCfg d=density_cfg(pipe); auto root=fs::absolute(fs::path("framework/results")/(cas+"_"+std::to_string(std::chrono::system_clock::now().time_since_epoch().count()))); fs::create_directories(root);
  std::ofstream summary(root/"case_summary.csv"); summary<<"stage,module,hpwl_before,hpwl_after,overflow_percent_before,overflow_percent_after,iterations,runtime_seconds\n";
  int idx=0; for(const auto& s:pipe.at("stages")) { auto mod=s.at("module").get<std::string>();
    if(compliance=="challenge_nonsmooth"&&mod=="historical_dct_poisson") throw std::runtime_error("challenge_nonsmooth rejects historical surrogate module");
    auto cfg=read_json(o.pipeline.parent_path()/s.at("config").get<std::string>()); M before=audit(db,d); auto t=std::chrono::steady_clock::now(); int its=run_module(mod,cfg,db,d,o.threads); double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count(); M after=audit(db,d);
    auto stage=root/("stage_"+(idx<10?std::string("0"):std::string())+std::to_string(idx)+"_"+mod); write_stage(stage,db,after,mod,idx,sec); csv(summary,idx,mod,before,after,sec,its); ++idx;
  }
  M final=audit(db,d); write_pl(db,root/"selected.pl"); Json manifest={{"case",cas},{"pipeline",esc(o.pipeline.string())},{"threads",o.threads},{"final_hpwl",final.hpwl},{"final_overflow_percent",final.overflow*100}}; std::ofstream(root/"run_manifest.json")<<manifest.dump(2)<<'\n'; std::cout<<"completed "<<cas<<" -> "<<root<<" HPWL="<<final.hpwl<<" overflow_percent="<<final.overflow*100<<"%"<<'\n';
}
}
int main(int ac,char**av){ if(ac>=2 && std::string(av[1])=="lab") return run_global_view_lab(ac,av); try{auto o=options(ac,av);if(o.command=="list-modules"){std::cout<<"layout_init\nhpwl_adam\nexact_joint_gp\nexact_recovery\nsurplus_bisection\nequal_shape_swap\nhistorical_dct_poisson (disabled)\nglobal_view_gp (lab)\n";return 0;}if(o.command=="audit"){if(o.placement.empty())throw std::runtime_error("audit requires --placement");ea::Database db=ea::read_bookshelf(o.dataset/o.cas/o.cas);ea::load_bookshelf_placement(db,o.placement);auto m=audit(db,{});std::cout<<std::setprecision(14)<<"hpwl="<<m.hpwl<<" density_energy="<<m.energy<<" overflow_percent="<<m.overflow*100<<"% max_density="<<m.maxd<<"\n";return 0;}if(o.command=="run"){run_case(o,o.cas);return 0;}if(o.command=="batch"&&o.all){for(auto&s:std::vector<std::string>{"adaptec1","adaptec2","adaptec3","adaptec4","bigblue1","bigblue2","bigblue3","bigblue4"})run_case(o,s);return 0;}usage();throw std::runtime_error("batch requires --all-cases");}catch(const std::exception&e){std::cerr<<"nsgp: "<<e.what()<<"\n";return 1;}}
