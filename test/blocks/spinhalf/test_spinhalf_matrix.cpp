#include "../../catch.hpp"

#include "testcases_spinhalf.h"
#include <hydra/all.h>
#include <iostream>

using namespace hydra;

TEST_CASE("spinhalf_matrix", "[spinhalf]") {
  using namespace hydra::testcases::spinhalf;

  {
    lila::Log.out("Spinhalf: Heisenberg chain test, J=1.0, N=2,..,6");
    for (int n_sites = 2; n_sites <= 6; ++n_sites)
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto [bonds, couplings, exact_eigs] =
            HBchain_fullspectrum_nup(n_sites, nup);
        auto block = Spinhalf<uint32_t>(n_sites, nup);
        auto H = MatrixReal(bonds, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));
        auto eigs = lila::EigenvaluesSym(H);
        REQUIRE(lila::close(eigs, exact_eigs));
      }
  }

  {
    lila::Log.out("Spinhalf: Heisenberg all-to-all tJ comparison");
    for (int n_sites = 2; n_sites <= 6; ++n_sites)
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto [bonds, couplings] = HB_alltoall(n_sites);

        auto block = Spinhalf<uint32_t>(n_sites, nup);
        auto block_tJ = tJ<uint32_t>(n_sites, nup, n_sites - nup);
        auto H = MatrixReal(bonds, couplings, block, block);
        auto H_tJ = MatrixReal(bonds, couplings, block_tJ, block_tJ);
        REQUIRE(lila::close(H, lila::Herm(H)));
        REQUIRE(lila::close(H_tJ, lila::Herm(H_tJ)));

        auto eigs = lila::EigenvaluesSym(H);
        auto eigs_tJ = lila::EigenvaluesSym(H_tJ);
        REQUIRE(lila::close(eigs, eigs_tJ));
      }
  }

  {
    lila::Log.out("Spinhalf: triangular N=12 complex exchange");
    int n_sites = 12;
    int nup = 6;
    std::vector<double> etas = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05};
    for (auto eta : etas) {
      auto [bonds, couplings, e0] = triangular_12_complex(nup, eta);
      auto block = Spinhalf<uint32_t>(n_sites, nup);
      auto H = MatrixCplx(bonds, couplings, block, block);
      REQUIRE(lila::close(H, lila::Herm(H)));

      auto eigs = lila::EigenvaluesSym(H);

      // comment: reference data from Lanczos, only ~10 digits precise
      lila::Log("eigs(0): {}, e0: {}", eigs(0), e0);
      REQUIRE(std::abs(eigs(0) - e0) < 1e-8);
    }
  }
}
