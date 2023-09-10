#include "../catch.hpp"

#include <hydra/combinatorics/lin_table.h>
#include <iostream>

template <class bit_t> void test_lintable(int n, int k) {
  using namespace hydra;
  combinatorics::LinTable<bit_t> lin_table(n, k);

  int64_t idx = 0;
  for (bit_t bits : combinatorics::Combinations<bit_t>(n, k)) {
    REQUIRE(lin_table.index(bits) == idx);
    ++idx;
  }
}

TEST_CASE("lintable", "[indexing]") {

  for (int n = 0; n < 10; ++n)
    for (int k = 0; k <= n; ++k) {
      test_lintable<uint16_t>(n, k);
      test_lintable<uint32_t>(n, k);
      test_lintable<uint64_t>(n, k);
    }
}
