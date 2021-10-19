#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

template <class bit_t>
void test_subsets(){
  using namespace hydra;
  using namespace hydra::combinatorics;
  using namespace hydra::utils;
  
  for (int n=0; n<8; ++n)
      {
	Subsets<bit_t> subs(n);
	REQUIRE(n == subs.n());
	
        idx_t ctr=0;
	bit_t current=0;
	for (auto sub : subs)
	  {
	    if (ctr != 0) REQUIRE(sub > current);
	    current = sub;
	    ++ctr;
	    REQUIRE(sub < ((bit_t)1 << n));
	  }
	REQUIRE(ctr == subs.size());
      }
}

TEST_CASE( "subsets", "[combinatorics/subsets]" ) {
  test_subsets<hydra::uint16>();
  test_subsets<hydra::uint32>();
  test_subsets<hydra::uint64>();
}
