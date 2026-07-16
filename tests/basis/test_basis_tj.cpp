// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/math/binomial.hpp>

using namespace xdiag;

// Sanity check of the (number-conserving) tJ basis enumeration: the size matches
// binomial(nsites,nup) * binomial(nsites-nup,ndn), and iterating yields exactly
// that many states, each satisfying the no-double-occupancy constraint with the
// right particle numbers (this also exercises the compressed -> full dn
// decompression in the iterator).
TEST_CASE("basis_tj", "[basis_tj]") {
  using enum_t = combinatorics::LinTable<uint32_t>;
  for (int64_t nsites = 1; nsites <= 6; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      for (int64_t ndn = 0; nup + ndn <= nsites; ++ndn) {
        enum_t enum_up(nsites, nup);
        enum_t enum_dncs(nsites - nup, ndn);
        basis::BasistJ<enum_t> b(enum_up, enum_dncs);

        int64_t expect =
            math::binomial(nsites, nup) * math::binomial(nsites - nup, ndn);
        REQUIRE(b.size() == expect);

        int64_t count = 0;
        for (auto [ups, dns] : b) {
          REQUIRE((ups & dns) == 0u);          // no double occupancy
          REQUIRE(bits::popcount(ups) == nup); // correct nup
          REQUIRE(bits::popcount(dns) == ndn); // correct ndn
          ++count;
        }
        REQUIRE(count == b.size());
      }
    }
  }
}
