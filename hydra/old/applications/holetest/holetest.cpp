#include <cstdlib>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>
#include <hydra/utils/combinatorics.h>
#include <lime/all.h>

lila::LoggerMPI lg;

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  // Checks down_hole_to_up implementation
  // Downspin config is 100101, "hole" configuration
  // is 10.
  // In this implementation, "hole" actually corresponds
  // to the number of remaining upspins, so the full 
  // configuration, in the form (upspins; downspins),
  // is (001000; 100101)
  
  uint64 downspins = 37;
  uint64 holes = 2;
  std::cout << "Upspin config is " << hydra::combinatorics::down_hole_to_up(downspins, holes) << std::endl;

  // Checks lin table indices against normal index table indices
  
  auto hole_hs_n4 = Spinhalf<uint64>(4, 2);
  auto lindexing_holes_n4 = LinTable<Spinhalf<uint64>, uint64>(hole_hs_n4);
  auto indexing_holes_n4 = IndexTable<Spinhalf<uint64>, uint64>(hole_hs_n4);
  
  for (auto state : hole_hs_n4) {
    std::cout << "Difference between normal index and lin table index is " 
      << lindexing_holes_n4.index(state) - indexing_holes_n4.index(state) << std::endl;
  }
  
auto hole_hs_n5 = Spinhalf<uint64>(5, 2);
  auto lindexing_holes_n5 = LinTable<Spinhalf<uint64>, uint64>(hole_hs_n5);
  auto indexing_holes_n5 = IndexTable<Spinhalf<uint64>, uint64>(hole_hs_n5);
  
  for (auto state : hole_hs_n5) {
    std::cout << "Difference between normal index and lin table index is " 
      << lindexing_holes_n5.index(state) - indexing_holes_n5.index(state) << std::endl;
  }
}
