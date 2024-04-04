#include "../catch.hpp"

#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/combinatorics/binomial.hpp>

TEST_CASE("bit_patterns", "[combinatorics]") {

  using namespace xdiag;
  using namespace xdiag::combinatorics;

  using bit_t = uint64_t;

  // Test spin pattern creation and indexing
  for (int n_sites = 1; n_sites < 6; ++n_sites)
    for (int n_upspins = 0; n_upspins <= n_sites; ++n_upspins) {
      bit_t state = ((bit_t)1 << n_upspins) - 1;
      for (int64_t n = 0; n < binomial(n_sites, n_upspins); ++n) {
        REQUIRE(state == get_nth_pattern<bit_t>(n, n_sites, n_upspins));
        REQUIRE(n == get_n_for_pattern(state, n_sites, n_upspins));
        state = get_next_pattern(state);
      }
    }
}
