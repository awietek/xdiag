#include "catch.hpp"

#include <iostream>
#include <hydra/all.h>

template <class state_t, class index_t>
void test_lintable_spinhalf(int n_sites, int n_upspins){
  using namespace hydra;
  using namespace hydra::hilbertspaces;
  using namespace hydra::indexing;

  using hs_t = Spinhalf<state_t>;
  hs_t hs(n_sites, n_upspins);
  LinTable<hs_t, index_t> lin_table(hs);
  IndexSearch<hs_t, index_t> index_search(hs);
  IndexSpinhalf<state_t, index_t> index_spinhalf(hs);

  index_t idx = 0;
  for (auto state : hs)
    {
      REQUIRE(lin_table.index(state) == idx);
      REQUIRE(lin_table.state(idx) == state); 
      REQUIRE(index_search.index(state) == idx);
      REQUIRE(index_search.state(idx) == state);
      REQUIRE(index_spinhalf.index(state) == idx);
      REQUIRE(index_spinhalf.state(idx) == state);
      ++idx;
    }
}

TEST_CASE( "indexing lintable test", "[indexing/lintable]" ) {
  using namespace hydra::all;

  for (int n_sites = 0; n_sites < 8; ++n_sites)
    for (int n_upspins = 0; n_upspins <= n_sites; ++n_upspins)
      {
	test_lintable_spinhalf<uint16, uint16>(n_sites, n_upspins);
	test_lintable_spinhalf<uint16, uint32>(n_sites, n_upspins);
	test_lintable_spinhalf<uint16, uint64>(n_sites, n_upspins);
	test_lintable_spinhalf<uint32, uint16>(n_sites, n_upspins);
	test_lintable_spinhalf<uint32, uint32>(n_sites, n_upspins);
	test_lintable_spinhalf<uint32, uint64>(n_sites, n_upspins);
	test_lintable_spinhalf<uint64, uint16>(n_sites, n_upspins);
	test_lintable_spinhalf<uint64, uint32>(n_sites, n_upspins);
	test_lintable_spinhalf<uint64, uint64>(n_sites, n_upspins);
      }
}
