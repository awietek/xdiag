// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/combinatorics/binomial.hpp>

TEST_CASE("bit_patterns", "[combinatorics]") {

  using namespace xdiag;
  using namespace xdiag::combinatorics;

  using bit_t = uint64_t;

  // Test spin pattern creation and indexing
  for (int nsites = 1; nsites < 6; ++nsites)
    for (int nupspins = 0; nupspins <= nsites; ++nupspins) {
      bit_t state = ((bit_t)1 << nupspins) - 1;
      for (int64_t n = 0; n < binomial(nsites, nupspins); ++n) {
        REQUIRE(state == get_nth_pattern<bit_t>(n, nsites, nupspins));
        REQUIRE(n == get_n_for_pattern(state, nsites, nupspins));
        state = get_next_pattern(state);
      }
    }
}
