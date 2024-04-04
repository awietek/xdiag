#include "../catch.hpp"

#include <iostream>
#include <xdiag/combinatorics/combinations.hpp>

template <typename bit_t> void test_combinations() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  for (int n = 0; n < 7; ++n) {
    for (int k = 0; k <= n; ++k) {
      Combinations<bit_t> combs(n, k);
      REQUIRE(n == combs.n());
      REQUIRE(k == combs.k());
      int64_t ctr = 0;
      bit_t current = 0;
      for (auto comb : combs) {

        if (ctr != 0)
          REQUIRE(comb > current);
        current = comb;
        ++ctr;
        REQUIRE(popcnt(comb) == k);
        REQUIRE(comb < ((bit_t)1 << n));
      }
      REQUIRE(ctr == combs.size());
    }
  }
}

TEST_CASE("Combinations", "[combinatorics]") {
  xdiag::Log.out("Testing Combinations");
  test_combinations<uint16_t>();
  test_combinations<uint32_t>();
  test_combinations<uint64_t>();
}
