#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

template <typename bit_t> void test_combinations_index() {
  using namespace hydra;
  using namespace hydra::combinatorics;
  using namespace hydra::bitops;

  for (int n = 0; n < 7; ++n)
    for (int k = 0; k <= n; ++k) {
      CombinationsIndex<bit_t> combs(n, k);
      REQUIRE(n == combs.n());
      REQUIRE(k == combs.k());
      idx_t ctr = 0;
      bit_t current = 0;
      for (auto [comb, idx] : combs) {

        if (ctr != 0)
          REQUIRE(comb > current);
        REQUIRE(ctr == idx);
        current = comb;
        ++ctr;
        REQUIRE(popcnt(comb) == k);
        REQUIRE(comb < ((bit_t)1 << n));
      }
      REQUIRE(ctr == combs.size());
    }

  // Test
  for (int n = 0; n < 7; ++n) {
    for (int k = 0; k <= n; ++k) {
      idx_t size = binomial(n, k);
      std::vector<bit_t> states_serial(size);
      std::vector<bit_t> states_parallel(size);

      // Create states serially
      for (auto [state, idx] : CombinationsIndex<bit_t>(n, k)) {
        states_serial[idx++] = state;
      }

      // Create states in parallel
#pragma openmp parallel
      {
        for (auto [state, idx] : CombinationsIndexThread<bit_t>(n, k)) {
          states_parallel[idx++] = state;
        }
      }

      for (idx_t idx = 0; idx < size; ++idx) {
        // Log("{} {}", states_serial[idx], states_parallel[idx]);
        REQUIRE(states_serial[idx] == states_parallel[idx]);
      }
    }
  }
}

TEST_CASE("CombinationsIndex", "[combinatorics]") {
  hydra::Log.out("Testing CombinationsIndex");
  test_combinations_index<uint16_t>();
  test_combinations_index<uint32_t>();
  test_combinations_index<uint64_t>();
}
