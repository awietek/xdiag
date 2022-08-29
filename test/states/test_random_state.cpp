#include "../catch.hpp"

#include "../blocks/electron/testcases_electron.h"
#include <hydra/all.h>
#include <iostream>
#include <set>

using namespace hydra;

TEST_CASE("random_state", "[states]") {
  using namespace hydra::testcases::electron;

  // Log.out("random state HB chain test");
  // for (int n_sites = 6; n_sites <= 8; ++n_sites) {
  //   Log("N={}", n_sites);
  //   auto [group, irreps] = get_cyclic_group_irreps(n_sites);

  //   for (int nup = 0; nup <= n_sites; ++nup) {
  //     for (auto irrep : irreps) {
  //       auto block = Spinhalf(n_sites, nup, group, irrep);
  //       auto state_real = RandomStateReal(block);

  //       auto set_real =
  //           std::set(state_real.vector().begin(), state_real.vector().end());
  //       REQUIRE(set_real.size() == state_real.size());
  //       LilaPrint(state_real.vector());

  //       auto state_cplx = RandomStateCplx(block);
  //       // auto set_cplx =
  //       //     std::set(state_cplx.vector().begin(),
  //       state_cplx.vector().end());
  //       // REQUIRE(set_cplx.size() == state_cplx.size());
  //       LilaPrint(state_cplx.vector());

  //     }
  //   }
  // }
#ifdef HYDRA_ENABLE_OPENMP

  // Check whether result with multiple threads is the same as with a single
  // thread
  auto block = Spinhalf(4);
  for (int seed = 0; seed < 10; ++seed) {
    auto state = RandomStateReal(block);
    auto state_cplx = RandomStateCplx(block);

    omp_set_num_threads(1);
    auto state2 = RandomStateReal(block);
    auto state_cplx2 = RandomStateCplx(block);
    REQUIRE(lila::Norm(state.vector() - state2.vector()) < 1e-12);
    REQUIRE(lila::Norm(state_cplx.vector() - state_cplx2.vector()) < 1e-12);
  }
#endif
}
