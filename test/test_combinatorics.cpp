#include "catch.hpp"

#include <hydra/all.h>

TEST_CASE( "utils/combinatorics test", "[utils/combinatorics]" ) {
  using uint64 = hydra::uint64;
  using hydra::combinatorics::binomial;
  using hydra::combinatorics::get_next_pattern;
  using hydra::combinatorics::get_nth_pattern;
  using hydra::combinatorics::get_n_for_pattern;

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

  // Test spin pattern creation and indexing
  for (int n_sites = 1; n_sites < 6; ++n_sites)
    for (int n_upspins = 0; n_upspins <= n_sites; ++n_upspins)
      {
	uint64 state = ((uint64)1 << n_upspins) - 1;
	for (int n = 0; n<binomial(n_sites, n_upspins); ++n)
	  {
	    // using hydra::hilbertspaces::PrintSpinhalf;
	    // std::cout << PrintSpinhalf(n_sites, state) << " "
	    // 	      << PrintSpinhalf(n_sites, get_nth_pattern((uint64)n, n_sites, n_upspins)) << " " 
	    // 	      << n << " " << get_n_for_pattern(state, n_sites, n_upspins) << std::endl;
	    REQUIRE(state == get_nth_pattern((uint64)n, n_sites, n_upspins));
	    REQUIRE(n == get_n_for_pattern(state, n_sites, n_upspins));
	    state = get_next_pattern(state);
	  }
      }

}
