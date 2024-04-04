#include "../catch.hpp"

template <class bit_t>
void test_up_down_hole(int n_sites, int n_upspins, int n_downspins) {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  using basis_t = BasisSpinHalf<bit_t>;
  {
    basis_t basisup(n_sites, n_upspins);
    basis_t basisholes(n_sites - n_upspins, n_sites - n_upspins - n_downspins);

    for (auto ups : basisup)
      for (auto holes : basisholes) {

        auto downs = up_hole_to_down(ups.spins, holes.spins);
        auto holes2 = up_down_to_hole(ups.spins, downs);
        REQUIRE(holes2 == holes.spins);
      }
  }

  {
    basis_t basisdown(n_sites, n_downspins);
    basis_t basisholes(n_sites - n_downspins,
                       n_sites - n_upspins - n_downspins);

    for (auto downs : basisdown)
      for (auto holes : basisholes) {
        auto ups = down_hole_to_up(downs.spins, holes.spins);
        auto holes2 = down_up_to_hole(downs.spins, ups);
        REQUIRE(holes2 == holes.spins);
      }
  }
}

TEST_CASE("combinatorics/up_down_hole", "[combinatorics]") {
  using namespace xdiag;
  
  for (int n_sites = 0; n_sites < 8; ++n_sites)
    for (int n_upspins = 0; n_upspins <= n_sites; ++n_upspins)
      for (int n_downspins = 0; n_downspins <= n_sites - n_upspins;
           ++n_downspins) {
	test_up_down_hole<uint16>(n_sites, n_upspins, n_downspins);
	test_up_down_hole<uint32>(n_sites, n_upspins, n_downspins);
	test_up_down_hole<xdiag::uint64>(n_sites, n_upspins, n_downspins);
      }
}
