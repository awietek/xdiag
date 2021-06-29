#include "../../catch.hpp"

#include "testcases_spinhalf.h"
#include <hydra/all.h>
#include <iostream>

using namespace hydra;

TEST_CASE("spinhalf_matrix", "[spinhalf]") {
  using namespace hydra::testcases::spinhalf;

  {
    HydraLog.out("Spinhalf: Heisenberg chain test, J=1.0, N=2,..,6");
    for (int n_sites = 2; n_sites <= 6; ++n_sites)
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto [bonds, couplings, exact_eigs] =
            HBchain_fullspectrum_nup(n_sites, nup);
        auto block = Spinhalf<uint32>(n_sites, nup);
        auto H = matrix_real(bonds, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));
        auto eigs = lila::EigenvaluesSym(H);
        REQUIRE(lila::close(eigs, exact_eigs));
      }
  }

  {
    HydraLog.out("Spinhalf: Heisenberg all-to-all tJ comparison");
    for (int n_sites = 2; n_sites <= 6; ++n_sites)
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto [bonds, couplings] = HB_alltoall(n_sites);

        auto block = Spinhalf<uint32>(n_sites, nup);
	auto block_tJ = tJ<uint32>(n_sites, nup, n_sites - nup);
        auto H = matrix_real(bonds, couplings, block, block);
        auto H_tJ = matrix_real(bonds, couplings, block_tJ, block_tJ);
        REQUIRE(lila::close(H, lila::Herm(H)));
	REQUIRE(lila::close(H_tJ, lila::Herm(H_tJ)));

        auto eigs = lila::EigenvaluesSym(H);
	auto eigs_tJ = lila::EigenvaluesSym(H_tJ);
        REQUIRE(lila::close(eigs, eigs_tJ));
      }
  }
}
