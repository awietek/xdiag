#include "../../catch.hpp"

#include <iostream>

#include "testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;

TEST_CASE("electron_apply", "[models]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;
  Couplings couplings;

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    printf("HubbardModel apply: random all-to-all test (real)\n");
    printf("N=%d\n", n_sites);
    std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);
    couplings["U"] = 5.0;

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron<uint32>(n_sites, nup, ndn);
        auto H = matrix_real(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<double>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = e0_real(bondlist, couplings, block);
        // HydraLog.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        REQUIRE(lila::close(e0_mat, e0_app));
      }
  }

  // /////////////////
  // Test Fermion all to all, free fermions (cplx, up/dn different)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    printf("HubbardModel: free fermion random all-to-all test, (cplx)\n");
    printf("N=%d\n", n_sites);
    std::tie(bondlist, couplings) = freefermion_alltoall_complex_updn(n_sites);
    couplings["U"] = 5.0;

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

	// Create block and matrix for comparison
        auto block = Electron<uint32>(n_sites, nup, ndn);
        auto H = matrix_cplx(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<complex>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = e0_cplx(bondlist, couplings, block);
        // HydraLog.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        REQUIRE(lila::close(e0_mat, e0_app));
      }
  }
}
