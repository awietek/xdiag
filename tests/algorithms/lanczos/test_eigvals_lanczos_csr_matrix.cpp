// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"

#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algebra/sparse/csr_matrix.hpp>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>

using namespace xdiag;

template <typename idx_t, typename coeff_t>
void test_eigvals_lanczos_csr_matrix() {
  using namespace xdiag::testcases::electron;

  OpSum ops;
  ops["U"] = 5.0;

  int nsites = 6;
  int num_eigenvalue = 1;

  if (std::is_same<coeff_t, double>::value) {
    printf("  eigvals_lanczos real test ...\n");
    ops = freefermion_alltoall(nsites);

    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {

        // Compute exact evals
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops, block, block);
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);

        // Compute evals with Lanczos
        auto spmat = csr_matrix<idx_t, coeff_t>(ops, block);
        auto res = eigvals_lanczos(ops, block, num_eigenvalue);
        auto evals = res.eigenvalues;
        double e0 = eigval0(spmat, block);
        REQUIRE(isapprox(evals_mat(0), e0));

        for (int i = 0; i < num_eigenvalue; ++i) {
          // lila::Log("a: {}, b: {}", evals_mat(i), evals_tmat(i));
          REQUIRE(std::abs(evals_mat(i) - evals(i)) < 1e-7);
        }
      }
  } else {
    printf("  eigvals_lanczos cplx test ...\n");
    ops = freefermion_alltoall_complex_updn(nsites);

    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron(nsites, nup, ndn);
        auto H = matrixC(ops, block, block);
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);

        // Compute evals with Lanczos
        auto spmat = csr_matrix<idx_t, coeff_t>(ops, block);
        auto res = eigvals_lanczos(spmat, block, num_eigenvalue);
        auto evals = res.eigenvalues;
        double e0 = eigval0(spmat, block);
        REQUIRE(isapprox(evals_mat(0), e0));

        for (int i = 0; i < num_eigenvalue; ++i)
          REQUIRE(isapprox(evals_mat(i), evals(i)));
      }
  }
}

TEST_CASE("eigvals_lanczos_csr_matrix", "[lanczos]") try {
  Log("eigvals_lanczos_csr_matrix int32_t double");
  test_eigvals_lanczos_csr_matrix<int32_t, double>();
  Log("eigvals_lanczos_csr_matrix int32_t complex");
  test_eigvals_lanczos_csr_matrix<int32_t, complex>();
  Log("eigvals_lanczos_csr_matrix int64_t double");
  test_eigvals_lanczos_csr_matrix<int64_t, double>();
  Log("eigvals_lanczos_csr_matrix int64_t complex");
  test_eigvals_lanczos_csr_matrix<int64_t, complex>();
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
