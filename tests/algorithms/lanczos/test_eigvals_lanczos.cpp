// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"

#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algebra/isapprox.hpp>

using namespace xdiag;

TEST_CASE("eigvals_lanczos", "[lanczos]") {
  using namespace xdiag::testcases::electron;

  OpSum ops;
  ops["U"] = 5.0;

  int nsites = 6;
  int num_eigenvalue = 1;

  printf("eigvals_lanczos real test ...\n");
  ops = freefermion_alltoall(nsites);

  for (int nup = 0; nup <= nsites; ++nup)
    for (int ndn = 0; ndn <= nsites; ++ndn) {

      // Compute exact evals
      auto block = Electron(nsites, nup, ndn);
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
  ops = freefermion_alltoall_complex_updn(nsites);

  for (int nup = 0; nup <= nsites; ++nup)
    for (int ndn = 0; ndn <= nsites; ++ndn) {

      // Create block and matrix for comparison
      auto block = Electron(nsites, nup, ndn);
      auto H = matrixC(ops, block, block);
      arma::vec evals_mat;
      arma::eig_sym(evals_mat, H);

      // Compute evals with Lanczos
      auto res = eigvals_lanczos(ops, block, num_eigenvalue);
      auto evals = res.eigenvalues;

      for (int i = 0; i < num_eigenvalue; ++i)
        REQUIRE(isapprox(evals_mat(i), evals(i)));
    }
  printf("Done.\n");
}
