#include "catch.hpp"

#include <iostream>

#include "combinatorics.h"
#include "hubbard.h"
#include "range.h"

template <class state_type>
void test_hubbard(){
  using hydra::combinatorics::binomial;
  using Hubbard = hydra::hilbertspaces::Hubbard<state_type>;
  using state_t = hydra::hilbertspaces::hubbard_state<state_type>;
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::hilbertspaces::Print;
  using hydra::hilbertspaces::Print;
  using hydra::utils::range;

  for (int n_sites : range<>(5))
    for (int n_upspins : range<>(n_sites+1))
      for (int n_downspins  : range<>(n_sites+1))
	{
	  state_t current = {0, 0};
	  hubbard_qn qn = {n_upspins, n_downspins};
	  // std::cout << n_sites << " " << qn.n_upspins << " " << qn.n_downspins << std::endl;
	  Hubbard hs(n_sites, qn);
	  int ctr=0;
	  for (auto state : hs)
	    {
	      // std::cout << Print(n_sites, state) << std::endl;
	      if (ctr!=0) REQUIRE(current < state);
	      current = state;
	      REQUIRE(current == state);
	      ++ctr;
	    }
	  // std::cout << std::endl;
	  REQUIRE(ctr == hs.size());
	}
}

TEST_CASE( "hilbertspaces hubbard test", "[hilbertspaces/hubbard]" ) {
  test_hubbard<hydra::uint16>();
  test_hubbard<hydra::uint32>();
  test_hubbard<hydra::uint64>();
}
