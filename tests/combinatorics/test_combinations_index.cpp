// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"
#include <xdiag/combinatorics/combinations_index.hpp>
#include <xdiag/utils/logger.hpp>

template <typename bit_t> void test_combinations_index() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  for (int n = 0; n < 7; ++n)
    for (int k = 0; k <= n; ++k) {
      CombinationsIndex<bit_t> combs(n, k);
      REQUIRE(n == combs.n());
      REQUIRE(k == combs.k());
      int64_t ctr = 0;
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
#ifdef _OPENMP
  for (int n = 0; n < 7; ++n) {
    for (int k = 0; k <= n; ++k) {
      int64_t size = binomial(n, k);
      std::vector<bit_t> states_serial(size);
      std::vector<bit_t> states_parallel(size);

      // Create states serially
      for (auto [state, idx] : CombinationsIndex<bit_t>(n, k)) {
        states_serial[idx++] = state;
      }

      // Create states in parallel
#pragma omp parallel
      {
        for (auto [state, idx] : CombinationsIndexThread<bit_t>(n, k)) {
          states_parallel[idx++] = state;
        }
      }

      for (int64_t idx = 0; idx < size; ++idx) {
        // Log("{} {}", states_serial[idx], states_parallel[idx]);
        REQUIRE(states_serial[idx] == states_parallel[idx]);
      }
    }
  }
#endif
}

TEST_CASE("CombinationsIndex", "[combinatorics]") {
  xdiag::Log("Testing CombinationsIndex");
  test_combinations_index<uint16_t>();
  test_combinations_index<uint32_t>();
  test_combinations_index<uint64_t>();
}
