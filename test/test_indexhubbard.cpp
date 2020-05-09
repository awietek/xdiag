#include "catch.hpp"

#include <iostream>

#include <hydra/all.h>


template <class indexing_t>
void test_indexhubbard(int n_sites, hydra::hilbertspaces::hubbard_qn qn){
  using namespace hydra;
  using namespace hydra::hilbertspaces;
  using namespace hydra::indexing;

  using index_t = typename indexing_t::index_t;
  using hs_t = Hubbard<typename indexing_t::state_t>;

  hs_t hs(n_sites, qn);
  IndexHubbard<indexing_t> index_hubbard(hs);

  index_t idx = 0;
  for (auto state : hs)
    {
      // std::cout << index_hubbard.index(state) << " " 
      // 		<< Print(n_sites, index_hubbard.state(idx)) << " " 
      // 		<< idx << " " 
      // 		<< Print(n_sites, state) << std::endl;
      REQUIRE(index_hubbard.index(state) == idx);
      REQUIRE(index_hubbard.state(idx) == state); 
      ++idx;
    }
  REQUIRE(idx == index_hubbard.size());
}

TEST_CASE( "indexing indexhubbard test", "[indexing/indexhubbard]" ) {
  using namespace hydra::all;

  for (int n_sites : range<>(5))
    for (int n_upspins : range<>(n_sites+1))
      for (int n_downspins  : range<>(n_sites+1))
      {
	hubbard_qn qn = {n_upspins, n_downspins};
	
	test_indexhubbard<IndexTable<Spinhalf<uint16>, uint16>>(n_sites, qn);
	test_indexhubbard<IndexTable<Spinhalf<uint32>, uint32>>(n_sites, qn);
	test_indexhubbard<IndexTable<Spinhalf<uint64>, uint64>>(n_sites, qn);

	test_indexhubbard<IndexSearch<Spinhalf<uint16>, uint16>>(n_sites, qn);
	test_indexhubbard<IndexSearch<Spinhalf<uint32>, uint32>>(n_sites, qn);
	test_indexhubbard<IndexSearch<Spinhalf<uint64>, uint64>>(n_sites, qn);

	test_indexhubbard<IndexSpinhalf<uint16, uint16>>(n_sites, qn);
	test_indexhubbard<IndexSpinhalf<uint32, uint32>>(n_sites, qn);
	test_indexhubbard<IndexSpinhalf<uint64, uint64>>(n_sites, qn);
 
      }
}
