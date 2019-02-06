#include "catch.hpp"

#include <iostream>

#include "charactertable.h"

TEST_CASE( "charactertable test", "[symmetries/charactertable]" ) {
  using namespace hydra::symmetries;

  auto character_table = read_charactertable("misc/lattice-files/square.16.spinlessfermions.pbc");
  Print(character_table);
}
