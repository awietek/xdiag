#include "../catch.hpp"

#include <hydra/all.h>

template <class bit_t> void test_electron() {
  using namespace hydra;

  for (int n_sites = 0; n_sites < 5; ++n_sites)
    for (int n_up = 0; n_up <= n_sites; ++n_up)
      for (int n_dn = 0; n_dn <= n_sites; ++n_dn) {
        state_electron<bit_t> current = {0, 0};
        qn_electron qn = {n_up, n_dn};

        BasisElectron<bit_t> basis(n_sites, qn);
        int ctr = 0;
        for (auto state : basis) {
          if (ctr != 0)
            REQUIRE(current < state);
          current = state;
          REQUIRE(current == state);
	  REQUIRE(valid(QN(state), n_sites));
          ++ctr;
        }
        REQUIRE(ctr == (int)basis.size());
      }
}

TEST_CASE("basis electron test", "[bases]") {
  test_electron<hydra::uint16>();
  test_electron<hydra::uint32>();
  test_electron<hydra::uint64>();
}
