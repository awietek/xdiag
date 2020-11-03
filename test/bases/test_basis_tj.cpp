#include "../catch.hpp"

#include <hydra/all.h>

template <class bit_t> void test_tj() {
  using namespace hydra;

  for (int n_sites = 0; n_sites < 7; ++n_sites)
    for (int n_up = 0; n_up <= n_sites; ++n_up)
      for (int n_dn = 0; n_dn <= n_sites - n_up; ++n_dn) {
        state_tj<bit_t> current = {0, 0};	
        BasisTJ<bit_t> basis(n_sites, {n_up, n_dn});
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

TEST_CASE("basis tj test", "[bases]") {
  test_tj<hydra::uint16>();
  test_tj<hydra::uint32>();
  test_tj<hydra::uint64>();
}
