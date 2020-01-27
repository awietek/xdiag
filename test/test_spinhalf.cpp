#include "catch.hpp"

#include <iostream>

#include <hydra/all.h>


template <class state_t>
void test_spinhalf(){
  using hydra::combinatorics::binomial;
  using Spinhalf = hydra::hilbertspaces::Spinhalf<state_t>;
  using hydra::hilbertspaces::PrintSpinhalf;
  using hydra::utils::range;

  for (int n_sites : range<>(6))
    for (int qn : range<>(n_sites+1))
      {
	state_t current = 0;
	Spinhalf hs(n_sites, qn);
	int ctr=0;
	for (auto state : hs)
	  {
	    if (ctr != 0) REQUIRE(state > current);
	    current = state;
	    // std::cout << qn << " " << PrintSpinhalf(n_sites, state) 
	    // 	      <<  " " << *Spinhalf(n_sites, qn).end() << std::endl;
	    ++ctr;
	  }
	REQUIRE(ctr == (int)hs.size());
      }
}

TEST_CASE( "hilbertspaces spinhalf test", "[hilbertspaces/spinhalf]" ) {
  test_spinhalf<hydra::uint16>();
  test_spinhalf<hydra::uint32>();
  test_spinhalf<hydra::uint64>();
}
