#include "../catch.hpp"

#include <hydra/all.h>

TEST_CASE( "combinatorics/bit_patterns", "[combinatorics]" ) {
  using namespace hydra;
  using namespace hydra::combinatorics;
  
  // Test spin pattern creation and indexing
  for (int n_sites = 1; n_sites < 6; ++n_sites)
    for (int n_upspins = 0; n_upspins <= n_sites; ++n_upspins)
      {
	uint64 state = ((uint64)1 << n_upspins) - 1;
	for (int n = 0; n<binomial(n_sites, n_upspins); ++n)
	  {
	    REQUIRE(state == get_nth_pattern((uint64)n, n_sites, n_upspins));
	    REQUIRE(n == get_n_for_pattern(state, n_sites, n_upspins));
	    state = get_next_pattern(state);
	  }
      }

}
