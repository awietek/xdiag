#include "catch.hpp"

#include <iostream>

#include <hydra/all.h>



TEST_CASE( "spacegroup test", "[symmetries/spacegroup]" ) {
  using namespace hydra::symmetries;

  std::vector<std::vector<int>> symmetries;
  symmetries.push_back({0, 1, 2, 3});
  symmetries.push_back({1, 2, 3, 0});
  symmetries.push_back({2, 3, 0, 1});
  symmetries.push_back({3, 0, 1, 2});
  
  SpaceGroup sg(symmetries);
  // HydraPrint(sg);
  
  SpaceGroup sub_sg = sg.subgroup({0, 2});
  // HydraPrint(sub_sg);
  
  // SpaceGroup square_sg = read_spacegroup("misc/lattice-files/square.16.spinlessfermions.pbc");
  // Print(square_sg);

}
