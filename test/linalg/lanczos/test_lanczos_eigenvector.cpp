#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;
using namespace lila;

TEST_CASE("lanczos_eigenvector", "[lanczos]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;
  Couplings couplings;
  couplings["U"] = 5.0;

  int n_sites = 6;
  int max_num_eigenvalue = 2;

  printf("LanczosEigenvectorReal test ...\n");
  std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);

  for (int nup = 2; nup <= n_sites / 2; ++nup)
    for (int ndn = 2; ndn <= n_sites / 2; ++ndn) {

      // Compute exact evals
      auto block = Electron(n_sites, nup, ndn);
      auto H = MatrixReal(bondlist, couplings, block, block);
      auto evals_mat = lila::EigenvaluesSym(H);

      // Compute evec with Lanczos
      for (int num_eigenvalue = 0; num_eigenvalue < max_num_eigenvalue;
           ++num_eigenvalue) {

        auto [tmat, evec] =
            LanczosEigenvectorReal(bondlist, couplings, block, num_eigenvalue);

        // Compute energy of eigenvector
        auto &v = evec;
        auto Hv = v;
        Apply(bondlist, couplings, block, v, block, Hv);
        auto e = Dot(v, Hv);
        REQUIRE(close(real(e), evals_mat(num_eigenvalue)));
        REQUIRE(close(imag(e), 0.0));
        REQUIRE(close(Norm(v), 1.0));
      }
    }
  printf("Done.\n");

  printf("LanczosEigenvectorCplx test ...\n");
  std::tie(bondlist, couplings) = freefermion_alltoall_complex_updn(n_sites);

  for (int nup = 2; nup <= n_sites / 3; ++nup)
    for (int ndn = 2; ndn <= n_sites / 3; ++ndn) {

      // Create block and matrix for comparison
      auto block = Electron(n_sites, nup, ndn);
      auto H = MatrixCplx(bondlist, couplings, block, block);
      auto evals_mat = lila::EigenvaluesSym(H);

      // Compute evec with Lanczos
      for (int num_eigenvalue = 0; num_eigenvalue < max_num_eigenvalue;
           ++num_eigenvalue) {

        auto [tmat, evec] =
            LanczosEigenvectorCplx(bondlist, couplings, block, num_eigenvalue);

        // Compute energy of eigenvector
        auto &v = evec;
        auto Hv = v;
        Apply(bondlist, couplings, block, v, block, Hv);
        auto e = Dot(v, Hv);
        REQUIRE(close(real(e), evals_mat(num_eigenvalue)));
        REQUIRE(close(imag(e), 0.0));
        REQUIRE(close(Norm(v), 1.0));
      }
    }
  printf("Done.\n");
}
