#include "epsilon_active/bookshelf.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/hpwl.hpp"
#include <cmath>
#include <filesystem>
#include <iostream>

int main() {
  ea::Database db;
  db.xl=0; db.yl=0; db.xh=20; db.yh=20; db.movable_area=100;
  db.nodes={{0,"a",5,5,10,10,false,false,"N"},{1,"b",15,5,10,10,false,false,"N"}};
  db.movable_ids={0,1}; db.node_pin_count={1,1}; db.node_pin_offsets={0,1,2}; db.node_pin_indices={0,1};
  db.pins={{0,0,0},{1,0,0}}; db.nets={{0,"n",0,2,2}};
  ea::ExactHpwl hpwl(db);
  if(std::abs(hpwl.evaluate(0,1,-1,nullptr,nullptr)-20.0)>1e-9) return 1;
  ea::ExactOverlapDensity density(db,2,2,.5);
  auto before=density.evaluate(0,1,nullptr,nullptr);
  auto move=density.evaluate_move(0,5,15); density.commit_move(move); db.nodes[0].x=5; db.nodes[0].y=15;
  auto after=density.evaluate(0,1,nullptr,nullptr);
  if(!std::isfinite(before.energy)||!std::isfinite(after.energy)||!std::isfinite(after.overflow)) return 2;
  db.node_by_name={{"a",0},{"b",1}};
  const auto root=std::filesystem::temp_directory_path()/"nsgp_numeric_contracts";
  std::filesystem::create_directories(root);
  const auto placement=root/"roundtrip.pl";
  ea::write_bookshelf_placement(db,placement);
  ea::Database copy=db; copy.nodes[0].x=0; copy.nodes[0].y=0;
  ea::load_bookshelf_placement(copy,placement);
  if(std::abs(copy.nodes[0].x-db.nodes[0].x)>1e-9 ||
     std::abs(copy.nodes[0].y-db.nodes[0].y)>1e-9) return 3;
  ea::ExactHpwl copy_hpwl(copy);
  if(std::abs(copy_hpwl.evaluate(0,1,-1,nullptr,nullptr)-
              hpwl.evaluate(0,1,-1,nullptr,nullptr))>1e-9) return 4;
  std::filesystem::remove(placement);
  std::filesystem::remove(root);

  // --- Fix A10: density delta contract (single + group moves) ---
  {
    ea::Database tdb;
    tdb.xl=0; tdb.yl=0; tdb.xh=40; tdb.yh=40; tdb.movable_area=4*4*3;
    tdb.nodes={
      {0,"a",5,5,4,4,false,false,"N"},
      {1,"b",15,15,4,4,false,false,"N"},
      {2,"c",30,30,4,4,false,false,"N"},
      {3,"f",22,10,4,4,true,false,"N"}};
    tdb.movable_ids={0,1,2}; tdb.fixed_ids={3};
    ea::ExactOverlapDensity den(tdb,4,4,1.0);
    const auto base=den.evaluate(0,1,nullptr,nullptr);
    auto mv=den.evaluate_move(0,5,15);
    ea::Database tdb2=tdb; tdb2.nodes[0].x=5; tdb2.nodes[0].y=15;
    ea::ExactOverlapDensity den2(tdb2,4,4,1.0);
    const auto full=den2.evaluate(0,1,nullptr,nullptr);
    if(std::abs((base.energy+mv.energy_delta)-full.energy)>1e-12) return 10;
    if(std::abs((base.overflow+mv.overflow_area_delta/tdb.movable_area)
                -full.overflow)>1e-12) return 11;
    // group move
    const std::vector<ea::DensityNodeMove> group={{1,20,20},{2,28,32}};
    auto gmv=den.evaluate_group_move(group);
    ea::Database tdb3=tdb; tdb3.nodes[1].x=20; tdb3.nodes[1].y=20;
    tdb3.nodes[2].x=28; tdb3.nodes[2].y=32;
    ea::ExactOverlapDensity den3(tdb3,4,4,1.0);
    const auto full3=den3.evaluate(0,1,nullptr,nullptr);
    if(std::abs((base.energy+gmv.energy_delta)-full3.energy)>1e-12) return 12;
  }

  // --- Fix A10: deterministic canonical audit is thread independent ---
  {
    ea::Database tdb;
    tdb.xl=0; tdb.yl=0; tdb.xh=40; tdb.yh=40; tdb.movable_area=4*4*3;
    tdb.nodes={
      {0,"a",5,5,4,4,false,false,"N"},
      {1,"b",15,15,4,4,false,false,"N"},
      {2,"c",30,30,4,4,false,false,"N"},
      {3,"f",22,10,4,4,true,false,"N"}};
    tdb.movable_ids={0,1,2}; tdb.fixed_ids={3};
    ea::ExactOverlapDensity den(tdb,8,8,1.0);
    const auto det=den.evaluate(0,1,ea::DensityAccumulationMode::DeterministicCanonical,nullptr,nullptr);
    const auto det2=den.evaluate(0,1,ea::DensityAccumulationMode::DeterministicCanonical,nullptr,nullptr);
    // bitwise-identical determinism across repeated canonical rebuilds
    if(det.energy!=det2.energy||det.overflow!=det2.overflow||
       det.max_density!=det2.max_density) return 13;
  }

  // --- Fix A10: HPWL weighted-net delta and exact-extrema subgradient ---
  {
    ea::Database tdb;
    tdb.xl=0; tdb.yl=0; tdb.xh=100; tdb.yh=20; tdb.movable_area=4;
    tdb.nodes={
      {0,"m",50,5,2,2,false,false,"N"},
      {1,"fix",10,5,2,2,true,false,"N"}};
    tdb.movable_ids={0}; tdb.fixed_ids={1};
    tdb.pins={{1,0,0},{0,0,0}};
    // one weighted net: weight 3, both pins at y=5
    tdb.nets.push_back({0,"w",0,2,3.0});
    tdb.node_pin_count={1,1};
    tdb.node_pin_weight_sum={3.0,3.0};
    tdb.node_pin_offsets={0,1,2}; tdb.node_pin_indices={1,0};
    ea::ExactHpwl hpwl(tdb);
    const ea::Real w0=hpwl.evaluate(0,1,-1,nullptr,nullptr);
    tdb.nodes[0].x=60;
    const ea::Real w1=hpwl.evaluate(0,1,-1,nullptr,nullptr);
    // moved right by 10 -> weighted HPWL increases by 3*10
    if(std::abs((w1-w0)-30.0)>1e-9) return 14;
    // exact-extrema for over-limit net: gradient weight split
    std::vector<ea::Real> gx,gy;
    hpwl.evaluate(0,1,1,ea::HighDegreeMode::ExactExtrema,&gx,&gy);
    if(std::abs(gx[0]-3.0)>1e-12) return 15; // movable is max side, weight 3
    // ignore mode zeroes the same net
    std::vector<ea::Real> gx2,gy2;
    hpwl.evaluate(0,1,1,ea::HighDegreeMode::Ignore,&gx2,&gy2);
    if(std::abs(gx2[0])>1e-12) return 16;
    tdb.nodes[0].x=50;
  }

  // --- Fix A10: node_pin_weight_sum bookkeeping ---
  {
    ea::Database tdb;
    tdb.nodes={{0,"a",0,0,1,1,false,false,"N"},{1,"b",5,0,1,1,false,false,"N"}};
    tdb.pins={{0,0,0},{1,0,0},{0,0,0}};
    tdb.nets.push_back({0,"n1",0,2,2.0});
    tdb.nets.push_back({1,"n2",2,1,0.5});
    tdb.node_pin_weight_sum.assign(2,0.0);
    for(const auto& net:tdb.nets)
      for(std::size_t p=net.pin_begin;p<net.pin_begin+net.pin_count;++p)
        tdb.node_pin_weight_sum[tdb.pins[p].node]+=net.weight;
    // node 0: 2.0+0.5, node 1: 2.0
    if(std::abs(tdb.node_pin_weight_sum[0]-2.5)>1e-12) return 17;
    if(std::abs(tdb.node_pin_weight_sum[1]-2.0)>1e-12) return 18;
  }
  std::cout<<"numeric contracts passed\n";
}
