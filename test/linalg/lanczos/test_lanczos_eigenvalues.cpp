#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;
using namespace lila;

TEST_CASE("lanczos_eigenvalues", "[lanczos]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;
  Couplings couplings;
  couplings["U"] = 5.0;

  int n_sites = 6;
  int num_eigenvalue = 2;

  printf("LanczosEigenvaluesReal test ...\n");
  std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);

  for (int nup = 2; nup <= n_sites / 2; ++nup)
    for (int ndn = 2; ndn <= n_sites / 2; ++ndn) {

      // Compute exact evals
      auto block = Electron<uint32>(n_sites, nup, ndn);
      auto H = MatrixReal(bondlist, couplings, block, block);
      auto evals_mat = lila::EigenvaluesSym(H);

      // Compute evals with Lanczos
      auto tmat =
          LanczosEigenvaluesReal(bondlist, couplings, block, num_eigenvalue);
      auto evals_tmat = tmat.eigenvalues();
      for (int i = 0; i < num_eigenvalue; ++i) {
	REQUIRE(close(evals_mat(i), evals_tmat(i)));
      }
    }
  printf("Done.\n");

  printf("LanczosEigenvaluesCplx test ...\n");
  std::tie(bondlist, couplings) = freefermion_alltoall_complex_updn(n_sites);

  for (int nup = 2; nup <= n_sites / 3; ++nup)
    for (int ndn = 2; ndn <= n_sites / 3; ++ndn) {

      // Create block and matrix for comparison
      auto block = Electron<uint32>(n_sites, nup, ndn);
      auto H = MatrixCplx(bondlist, couplings, block, block);
      auto evals_mat = lila::EigenvaluesSym(H);

      // Compute evals with Lanczos
      auto tmat =
          LanczosEigenvaluesCplx(bondlist, couplings, block, num_eigenvalue);
      auto evals_tmat = tmat.eigenvalues();
      for (int i = 0; i < num_eigenvalue; ++i)
        REQUIRE(close(evals_mat(i), evals_tmat(i)));
    }
  printf("Done.\n");

}
