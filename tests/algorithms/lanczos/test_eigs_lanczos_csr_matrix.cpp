// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algebra/sparse/csr_matrix.hpp>
#include <xdiag/algorithms/lanczos/eigs_lanczos.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>

#include <xdiag/algebra/isapprox.hpp>

using namespace xdiag;

template <typename idx_t, typename coeff_t>
void test_eigs_lanczos_csr_matrix() {
  using namespace xdiag::testcases::electron;

  OpSum ops;

  int nsites = 6;
  int max_num_eigenvalue = 2;

  if (std::is_same<coeff_t, double>::value) {
    printf("  lanczos_eigenvector_real test ...\n");
    ops = freefermion_alltoall(nsites);
    ops["U"] = 5.0;

    for (int nup = 2; nup <= nsites / 2; ++nup)
      for (int ndn = 2; ndn <= nsites / 2; ++ndn) {

        // Compute exact evals
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops, block, block);
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);

        // Compute evec with Lanczos
        auto spmat = csr_matrix<idx_t, coeff_t>(ops, block);
        auto res = eigs_lanczos(spmat, block, max_num_eigenvalue);
        auto evals = res.eigenvalues;
        auto [e0, _] = eig0(spmat, block);
        REQUIRE(isapprox(evals_mat(0), e0));

        // Compute energy of eigenvector
        for (int num_eigenvalue = 0; num_eigenvalue < max_num_eigenvalue;
             ++num_eigenvalue) {
          auto v = res.eigenvectors.col(num_eigenvalue);
          // XDIAG_SHOW(norm(v));
          REQUIRE(isapprox(norm(v), 1.0));

          auto Hv = v;
          apply(ops, v, Hv);
          double e = dot(v, Hv);

          // XDIAG_SHOW(num_eigenvalue);
          // XDIAG_SHOW(e);
          // XDIAG_SHOW(evals_mat(num_eigenvalue));
          REQUIRE(isapprox(real(e), evals_mat(num_eigenvalue)));
          REQUIRE(isapprox(imag(e), 0.0));
        }
      }
  }

  if (std::is_same<coeff_t, complex>::value) {

    printf("  lanczos_eigenvector_cplx test ...\n");
    ops = freefermion_alltoall_complex_updn(nsites);
    for (int nup = 2; nup <= nsites / 3; ++nup)
      for (int ndn = 2; ndn <= nsites / 3; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron(nsites, nup, ndn);
        auto H = matrixC(ops, block, block);
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);

        // Compute evec with Lanczos

        auto spmat = csr_matrix<idx_t, coeff_t>(ops, block);
        auto res =
            eigs_lanczos(spmat, block, max_num_eigenvalue, 1e-12, 1000, true);
        auto [e0, _] = eig0(spmat, block);
        REQUIRE(isapprox(evals_mat(0), e0));

        // Compute energy of eigenvector
        for (int num_eigenvalue = 0; num_eigenvalue < max_num_eigenvalue;
             ++num_eigenvalue) {
          auto v = res.eigenvectors.col(num_eigenvalue);
          REQUIRE(isapprox(norm(v), 1.0));
          auto Hv = v;
          apply(ops, v, Hv);
          auto e = dotC(v, Hv);
          REQUIRE(isapprox(real(e), evals_mat(num_eigenvalue)));
          REQUIRE(isapprox(imag(e), 0.0));
        }
      }
  }
}

TEST_CASE("eigs_lanczos_csr_matrix", "[lanczos]") try {
  Log("eigs_lanczos_csr_matrix int32_t double");
  test_eigs_lanczos_csr_matrix<int32_t, double>();
  Log("eigs_lanczos_csr_matrix int32_t complex");
  test_eigs_lanczos_csr_matrix<int32_t, complex>();
  Log("eigs_lanczos_csr_matrix int64_t double");
  test_eigs_lanczos_csr_matrix<int64_t, double>();
  Log("eigs_lanczos_csr_matrix int64_t complex");
  test_eigs_lanczos_csr_matrix<int64_t, complex>();
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
