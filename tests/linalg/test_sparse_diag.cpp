// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/blocks/spinhalf/testcases_spinhalf.hpp>
#include <tests/catch.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/linalg/sparse_diag.hpp>

using namespace xdiag;

// Smoke test that every entry point in sparse_diag.hpp is wired up correctly
// (eigval0/eig0 -> Lanczos, eigvals/eigs -> LOBPCG), for both OpSum and
// CSRMatrix inputs. Small block, no heavy computation.
TEST_CASE("sparse_diag", "[linalg]") {
  using namespace xdiag::testcases::spinhalf;

  int nsites = 8;
  int nup = 4;
  OpSum ops = HBchain(nsites, 1.0);
  auto block = Spinhalf(nsites, nup);

  // Dense reference.
  arma::mat H = matrix(ops, block);
  arma::vec evals_ref;
  arma::eig_sym(evals_ref, H);

  int64_t neigs = 2;

  auto residual_ok = [&](State const &v, arma::vec const &vals) {
    arma::mat V = v.matrix();
    arma::mat AV(V.n_rows, V.n_cols);
    apply(ops, block, V, block, AV);
    for (arma::uword i = 0; i < V.n_cols; ++i) {
      REQUIRE(arma::norm(AV.col(i) - vals(i) * V.col(i)) < 1e-6);
    }
  };

  // --- OpSum ---
  SECTION("OpSum") {
    // eigval0 -> Lanczos ground energy
    double e0 = eigval0(ops, block);
    REQUIRE(std::abs(e0 - evals_ref(0)) < 1e-8);

    // eig0 -> ground energy + eigenvector
    auto [e0b, psi0] = eig0(ops, block);
    REQUIRE(std::abs(e0b - evals_ref(0)) < 1e-8);
    residual_ok(psi0, arma::vec{e0b});

    // eigvals -> LOBPCG lowest neigs
    arma::vec vals = eigvals(ops, block, neigs);
    REQUIRE((int64_t)vals.size() == neigs);
    for (int64_t i = 0; i < neigs; ++i) {
      REQUIRE(std::abs(vals(i) - evals_ref(i)) < 1e-6);
    }

    // eigs -> LOBPCG lowest neigs + eigenvectors
    auto [vals2, psis] = eigs(ops, block, neigs);
    REQUIRE((int64_t)vals2.size() == neigs);
    REQUIRE(psis.ncols() == neigs);
    for (int64_t i = 0; i < neigs; ++i) {
      REQUIRE(std::abs(vals2(i) - evals_ref(i)) < 1e-6);
    }
    residual_ok(psis, vals2);
  }

  // --- CSRMatrix ---
  SECTION("CSRMatrix") {
    CSRMatrix<int64_t, double> A = csr_matrix(ops, block);
    A.ishermitian = true;

    REQUIRE(std::abs(eigval0(A, block) - evals_ref(0)) < 1e-8);

    auto [e0, psi0] = eig0(A, block);
    REQUIRE(std::abs(e0 - evals_ref(0)) < 1e-8);

    arma::vec vals = eigvals(A, block, neigs);
    REQUIRE((int64_t)vals.size() == neigs);
    for (int64_t i = 0; i < neigs; ++i) {
      REQUIRE(std::abs(vals(i) - evals_ref(i)) < 1e-6);
    }

    auto [vals2, psis] = eigs(A, block, neigs);
    REQUIRE(psis.ncols() == neigs);
    for (int64_t i = 0; i < neigs; ++i) {
      REQUIRE(std::abs(vals2(i) - evals_ref(i)) < 1e-6);
    }
  }
}
