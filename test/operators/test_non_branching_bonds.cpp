#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

TEST_CASE("non_branching_bonds", "[operators]") {
  using namespace hydra;

  arma::mat sx = {{0., 1.}, {1., 0.}};
  auto sx_bond = Bond(sx, 1);
  REQUIRE(operators::is_non_branching_bond(sx_bond));
}
