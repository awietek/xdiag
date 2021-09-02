#include "../catch.hpp"

#include <iostream>
#include <hydra/all.h>

template <class bit_t>
void test_lintable(int n, int k){
  using namespace hydra;
  
  indexing::LinTable<bit_t> lin_table(n, k);
  
  idx_t idx = 0;
  for (bit_t bits : Combinations<bit_t>(n, k))
    {
      REQUIRE(lin_table.index(bits) == idx);
      ++idx;
    }
}

TEST_CASE( "lintable", "[indexing]" ) {
  
  for (int n = 0; n < 10; ++n)
    for (int k = 0; k <= n; ++k)
      {
	test_lintable<hydra::uint16>(n, k);
	test_lintable<hydra::uint32>(n, k);
	test_lintable<hydra::uint64>(n, k);
      }

}
