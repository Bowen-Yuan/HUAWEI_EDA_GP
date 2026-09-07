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
  std::cout<<"numeric contracts passed\n";
}
