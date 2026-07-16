// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <iostream>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/catch.hpp>

#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/linalg/lanczos/eigvals_lanczos.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("eigvals_lanczos", "[lanczos]") {
  using namespace xdiag::testcases::electron;

  OpSum ops;
  ops["U"] = 5.0;

  int nsites = 6;
  int num_eigenvalue = 1;

  Log("eigvals_lanczos real test ...");
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
          REQUIRE(std::abs(evals_mat(i) - evals(i)) < 1e-7);
        }
      } catch (Error e) {
        error_trace(e);
        throw;
      }
    }
  Log("Done.");

  Log("eigvals_lanczos cplx test ...");
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
  Log("Done.");
}
