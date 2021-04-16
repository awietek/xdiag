#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

template <class bit_t>
void test_combinations(){
  using namespace hydra;
  using namespace hydra::combinatorics;
  using namespace hydra::utils;
  
  for (int n=0; n<7; ++n)
    for (int k=0; k<=n; ++k)
      {
	Combinations<bit_t> combs(n, k);
	REQUIRE(n == combs.n());
	REQUIRE(k == combs.k());
	
        idx_t ctr=0;
	bit_t current=0;
	for (auto comb : combs)
	  {
	    if (ctr != 0) REQUIRE(comb > current);
	    current = comb;
	    ++ctr;
	    REQUIRE(popcnt(comb) == k);
	    REQUIRE(comb < (1 << n));
	  }
	REQUIRE(ctr == combs.size());
      }
}

TEST_CASE( "combinations", "[combinatorics/combinations]" ) {
  test_combinations<hydra::uint16>();
  test_combinations<hydra::uint32>();
  test_combinations<hydra::uint64>();
}
