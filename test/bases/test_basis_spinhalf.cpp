#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>


template <class bit_t>
void test_spinhalf(){
  using namespace hydra;
  
  for (int n_sites=2; n_sites<7; ++n_sites)
    for (int n_up=0; n_up<=n_sites; ++n_up)
      {
        bit_t current = 0;
	BasisSpinHalf<bit_t> basis(n_sites, n_up);
	int ctr=0;
	for (auto state : basis)
	  {
	    if (ctr != 0) REQUIRE(state.spins > current);
	    current = state.spins;
	    ++ctr;
	    REQUIRE(valid(QN(state), n_sites));
	  }
	REQUIRE(ctr == (int)basis.size());
      }
}

TEST_CASE( "hilbertspaces spinhalf test", "[hilbertspaces/spinhalf]" ) {
  test_spinhalf<hydra::uint16>();
  test_spinhalf<hydra::uint32>();
  test_spinhalf<hydra::uint64>();
}
