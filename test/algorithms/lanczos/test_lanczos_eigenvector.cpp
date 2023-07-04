#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.h"
#include <hydra/algorithms/lanczos/lanczos_eigenvector.h>
#include <hydra/blocks/electron/electron_matrix.h>
#include <hydra/utils/close.h>

using namespace hydra;

TEST_CASE("lanczos_eigenvector", "[lanczos]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;

  int n_sites = 6;
  int max_num_eigenvalue = 2;

  printf("lanczos_eigenvector_real test ...\n");
  bondlist = freefermion_alltoall(n_sites);
  bondlist["U"] = 5.0;

  for (int nup = 2; nup <= n_sites / 2; ++nup)
    for (int ndn = 2; ndn <= n_sites / 2; ++ndn) {

      // Compute exact evals
      auto block = Electron(n_sites, nup, ndn);
      auto H = matrix_real(bondlist, block, block);
      arma::vec evals_mat;
      arma::eig_sym(evals_mat, H);

      // Compute evec with Lanczos
      for (int num_eigenvalue = 0; num_eigenvalue < max_num_eigenvalue;
           ++num_eigenvalue) {

        auto [tmat, evec] =
            lanczos_eigenvector_real(bondlist, block, num_eigenvalue);

        (void)tmat;

        // Compute energy of eigenvector
        auto &v = evec;
        auto Hv = v;
        apply(bondlist, block, v, block, Hv);
        auto e = arma::cdot(v, Hv);
        REQUIRE(close(real(e), evals_mat(num_eigenvalue)));
        REQUIRE(close(imag(e), 0.0));
        REQUIRE(close(arma::norm(v), 1.0));
      }
    }
  printf("Done.\n");

  printf("lanczos_eigenvector_cplx test ...\n");
  bondlist = freefermion_alltoall_complex_updn(n_sites);

  for (int nup = 2; nup <= n_sites / 3; ++nup)
    for (int ndn = 2; ndn <= n_sites / 3; ++ndn) {

      // Create block and matrix for comparison
      auto block = Electron(n_sites, nup, ndn);
      auto H = matrix_cplx(bondlist, block, block);
      arma::vec evals_mat;
      arma::eig_sym(evals_mat, H);

      // Compute evec with Lanczos
      for (int num_eigenvalue = 0; num_eigenvalue < max_num_eigenvalue;
           ++num_eigenvalue) {

        auto [tmat, evec] =
            lanczos_eigenvector_cplx(bondlist, block, num_eigenvalue);
        (void)tmat;
	
        // Compute energy of eigenvector
        auto &v = evec;
        auto Hv = v;
        apply(bondlist, block, v, block, Hv);
        auto e = arma::cdot(v, Hv);
        REQUIRE(close(real(e), evals_mat(num_eigenvalue)));
        REQUIRE(close(imag(e), 0.0));
        REQUIRE(close(arma::norm(v), 1.0));
      }
    }
  printf("Done.\n");
}
