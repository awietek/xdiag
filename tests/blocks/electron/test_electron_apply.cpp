// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <vector>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/blocks/tj/testcases_tj.hpp>
#include <tests/catch.hpp>
#include <tests/is_approx_hermitian.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// For a number-non-conserving Hamiltonian the full-block ground state must equal
// the minimum over all (nup, ndn) number sectors.
static void test_electron_np_no_np_apply(int nsites, OpSum ops) {
  auto block_full = Electron(nsites);
  auto e0_full = eigval0(ops, block_full);

  std::vector<double> e0s;
  for (int nup = 0; nup <= nsites; ++nup) {
    for (int ndn = 0; ndn <= nsites; ++ndn) {
      auto block = Electron(nsites, nup, ndn);
      e0s.push_back(eigval0(ops, block));
    }
  }
  auto e0_np = *std::min_element(e0s.begin(), e0s.end());
  REQUIRE(isapprox(e0_full, e0_np));
}

TEST_CASE("electron_apply", "[electron]") try {
  using namespace xdiag::testcases::electron;

  OpSum ops;

  // Hubbard random all-to-all, free fermions (real)
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_apply: Hubbard random all-to-all test (real), N: {}", nsites);
    ops = freefermion_alltoall(nsites);
    ops["U"] = 5.0;

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops, block, block);
        REQUIRE(testcases::is_approx_hermitian(H, 1e-8));

        // apply matches dense matrix multiplication (vector and multi-column)
        arma::cx_vec v(block.size(), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(block.size(), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));
        arma::mat m(block.size(), 5, arma::fill::randn);
        arma::mat n1 = H * m;
        arma::mat n2(block.size(), 5, arma::fill::zeros);
        apply(ops, block, m, block, n2);
        REQUIRE(isapprox(n1, n2));

        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_app = eigval0(ops, block);
        CHECK(isapprox(evals_mat(0), e0_app, 1e-8, 1e-8));
      }
    }
  }

  // Hubbard random all-to-all, free fermions (cplx, up/dn different)
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_apply: Hubbard random all-to-all test (cplx), N: {}", nsites);
    ops = freefermion_alltoall_complex_updn(nsites);
    ops["U"] = 5.0;

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        auto H = matrixC(ops, block, block);
        REQUIRE(testcases::is_approx_hermitian(H, 1e-8));

        arma::cx_vec v(block.size(), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(block.size(), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));

        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_app = eigval0(ops, block);
        CHECK(isapprox(evals_mat(0), e0_app, 1e-8, 1e-8));
      }
    }
  }

  // Henry's MATLAB code test (Hubbard + hopping + Heisenberg terms)
  Log("electron_apply: U-hopping-HB apply of Henry's Matlab code");
  {
    int nsites = 4;
    auto [ops4, eigs_correct] = randomAlltoAll4NoU();
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops4, block, block);
        REQUIRE(testcases::is_approx_hermitian(H, 1e-8));
        arma::vec v(block.size(), arma::fill::randn);
        arma::vec w1 = H * v;
        arma::vec w2(block.size(), arma::fill::zeros);
        apply(ops4, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        CHECK(isapprox(evals_mat(0), eigval0(ops4, block), 1e-8, 1e-8));
      }
    }

    auto [opsU, eigs_correctU] = randomAlltoAll4();
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(opsU, block, block);
        REQUIRE(testcases::is_approx_hermitian(H, 1e-8));
        arma::vec v(block.size(), arma::fill::randn);
        arma::vec w1 = H * v;
        arma::vec w2(block.size(), arma::fill::zeros);
        apply(opsU, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        CHECK(isapprox(evals_mat(0), eigval0(opsU, block), 1e-6, 1e-6));
      }
    }
  }

  // Number-conserving vs full block, complex tJ-style all-to-all
  for (int N = 3; N <= 5; ++N) {
    Log("electron_apply: random all-to-all complex exchange Np <-> NoNp, N={}",
        N);
    auto ops_tj = xdiag::testcases::tj::tj_alltoall_complex(N);
    test_electron_np_no_np_apply(N, ops_tj);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
