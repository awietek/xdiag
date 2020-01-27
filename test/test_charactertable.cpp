#include "catch.hpp"

#include <iostream>

#include <hydra/all.h>

TEST_CASE( "charactertable test", "[symmetries/charactertable]" ) {
  using namespace hydra::symmetries;

  auto character_table = read_charactertable("../misc/lattice-files/square.16.spinlessfermions.pbc");
  // HydraPrint(character_table);
}
