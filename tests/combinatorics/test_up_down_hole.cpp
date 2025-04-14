// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

template <class bit_t>
void test_up_down_hole(int nsites, int nupspins, int n_downspins) {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  using basis_t = BasisSpinHalf<bit_t>;
  {
    basis_t basisup(nsites, nupspins);
    basis_t basisholes(nsites - nupspins, nsites - nupspins - n_downspins);

    for (auto ups : basisup)
      for (auto holes : basisholes) {

        auto downs = up_hole_to_down(ups.spins, holes.spins);
        auto holes2 = up_down_to_hole(ups.spins, downs);
        REQUIRE(holes2 == holes.spins);
      }
  }

  {
    basis_t basisdown(nsites, n_downspins);
    basis_t basisholes(nsites - n_downspins,
                       nsites - nupspins - n_downspins);

    for (auto downs : basisdown)
      for (auto holes : basisholes) {
        auto ups = down_hole_to_up(downs.spins, holes.spins);
        auto holes2 = downup_to_hole(downs.spins, ups);
        REQUIRE(holes2 == holes.spins);
      }
  }
}

TEST_CASE("combinatorics/up_down_hole", "[combinatorics]") {
  using namespace xdiag;
  
  for (int nsites = 0; nsites < 8; ++nsites)
    for (int nupspins = 0; nupspins <= nsites; ++nupspins)
      for (int n_downspins = 0; n_downspins <= nsites - nupspins;
           ++n_downspins) {
	test_up_down_hole<uint16>(nsites, nupspins, n_downspins);
	test_up_down_hole<uint32>(nsites, nupspins, n_downspins);
	test_up_down_hole<xdiag::uint64>(nsites, nupspins, n_downspins);
      }
}
