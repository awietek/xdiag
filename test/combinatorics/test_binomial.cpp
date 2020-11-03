#include "../catch.hpp"

#include <hydra/all.h>

TEST_CASE( "combinatorics/binomial", "[combinatorics]" ) {
  using namespace hydra::combinatorics;

  REQUIRE(binomial(1, -1) == 0);
  REQUIRE(binomial(1, 0) == 1);
  REQUIRE(binomial(1, 1) == 1);
  REQUIRE(binomial(1, 2) == 0);

  REQUIRE(binomial(2, -1) == 0);
  REQUIRE(binomial(2, 0) == 1);
  REQUIRE(binomial(2, 1) == 2);
  REQUIRE(binomial(2, 2) == 1);
  REQUIRE(binomial(2, 3) == 0);

  REQUIRE(binomial(3, -1) == 0);
  REQUIRE(binomial(3, 0) == 1);
  REQUIRE(binomial(3, 1) == 3);
  REQUIRE(binomial(3, 2) == 3);
  REQUIRE(binomial(3, 3) == 1);
  REQUIRE(binomial(3, 4) == 0);

  REQUIRE(binomial(4, -1) == 0);
  REQUIRE(binomial(4, 0) == 1);
  REQUIRE(binomial(4, 1) == 4);
  REQUIRE(binomial(4, 2) == 6);
  REQUIRE(binomial(4, 3) == 4);
  REQUIRE(binomial(4, 4) == 1);
  REQUIRE(binomial(4, 5) == 0);

}
