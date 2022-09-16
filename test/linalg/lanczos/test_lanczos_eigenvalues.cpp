#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;

TEST_CASE("lanczos_eigenvalues", "[lanczos]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;
  Couplings couplings;
  couplings["U"] = 5.0;

  int n_sites = 6;
  int num_eigenvalue = 1;

  printf("LanczosEigenvaluesReal test ...\n");
  std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);

  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      // Compute exact evals
      auto block = Electron(n_sites, nup, ndn);
      auto H = MatrixReal(bondlist, couplings, block, block);
      arma::vec evals_mat;
      arma::eig_sym(evals_mat, H);

      // Compute evals with Lanczos
      auto tmat =
	LanczosEigenvaluesReal(bondlist, couplings, block, num_eigenvalue);
      auto evals_tmat = tmat.eigenvalues();
      for (int i = 0; i < num_eigenvalue; ++i) {
	// lila::Log("a: {}, b: {}", evals_mat(i), evals_tmat(i));
        REQUIRE(close(evals_mat(i), evals_tmat(i)));
      }
    }
  printf("Done.\n");

  printf("LanczosEigenvaluesCplx test ...\n");
  std::tie(bondlist, couplings) = freefermion_alltoall_complex_updn(n_sites);

  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      // Create block and matrix for comparison
      auto block = Electron(n_sites, nup, ndn);
      auto H = MatrixCplx(bondlist, couplings, block, block);
      arma::vec evals_mat;
      arma::eig_sym(evals_mat, H);
      
      // Compute evals with Lanczos
      auto tmat =
          LanczosEigenvaluesCplx(bondlist, couplings, block, num_eigenvalue);
      auto evals_tmat = tmat.eigenvalues();
      
      for (int i = 0; i < num_eigenvalue; ++i)
        REQUIRE(close(evals_mat(i), evals_tmat(i)));
    }
  printf("Done.\n");
}
