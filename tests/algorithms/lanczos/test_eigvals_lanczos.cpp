#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"

#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/blocks/electron/electron_matrix.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;

TEST_CASE("eigvals_lanczos", "[lanczos]") {
  using namespace xdiag::testcases::electron;

  OpSum ops;
  ops["U"] = 5.0;

  int n_sites = 6;
  int num_eigenvalue = 1;

  printf("eigvals_lanczos real test ...\n");
  ops = freefermion_alltoall(n_sites);

  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      // Compute exact evals
      auto block = Electron(n_sites, nup, ndn);
      auto H = matrix(ops, block, block);
      arma::vec evals_mat;
      arma::eig_sym(evals_mat, H);

      // Compute evals with Lanczos
      try {
        auto res = eigvals_lanczos(ops, block, num_eigenvalue);
        auto evals = res.eigenvalues;
        for (int i = 0; i < num_eigenvalue; ++i) {
          // lila::Log("a: {}, b: {}", evals_mat(i), evals_tmat(i));
          REQUIRE(std::abs(evals_mat(i) - evals(i)) < 1e-7 );
        }
      } catch (Error e) {
        error_trace(e);
      }
    }
  printf("Done.\n");

  printf("eigvals_lanczos cplx test ...\n");
  ops = freefermion_alltoall_complex_updn(n_sites);

  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      // Create block and matrix for comparison
      auto block = Electron(n_sites, nup, ndn);
      auto H = matrixC(ops, block, block);
      arma::vec evals_mat;
      arma::eig_sym(evals_mat, H);

      // Compute evals with Lanczos
      auto res = eigvals_lanczos(ops, block, num_eigenvalue);
      auto evals = res.eigenvalues;

      for (int i = 0; i < num_eigenvalue; ++i)
        REQUIRE(close(evals_mat(i), evals(i)));
    }
  printf("Done.\n");
}
