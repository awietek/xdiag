// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>
#include <tests/blocks/electron/testcases_electron.hpp>

#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <tests/blocks/spinhalf/testcases_spinhalf.hpp>
#include <xdiag/linalg/lobpcg/eigs_lobpcg.hpp>
#include <xdiag/linalg/lobpcg/lobpcg.hpp>
#include <xdiag/math/dot.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/xdiag_show.hpp>
#include <xdiag/utils/logger.hpp>

#include <cmath>
#include <vector>

using namespace xdiag;

// Groups sorted eigenvalues into degenerate clusters and returns the size of
// each cluster (the multiplicities).
static std::vector<int64_t> multiplicities(arma::vec const &e, double tol) {
  std::vector<int64_t> mult;
  int64_t i = 0;
  while (i < (int64_t)e.n_elem) {
    int64_t j = i + 1;
    while ((j < (int64_t)e.n_elem) && (std::abs(e(j) - e(i)) < tol)) {
      ++j;
    }
    mult.push_back(j - i);
    i = j;
  }
  return mult;
}

TEST_CASE("lobpcg", "[linalg]") {
  using namespace xdiag::testcases::electron;

  int nsites = 6;
  int nup = 3;
  int ndn = 3;
  OpSum ops = freefermion_alltoall(nsites);
  auto block = Electron(nsites, nup, ndn);

  // Dense reference spectrum.
  arma::mat H = matrix(ops, block);
  arma::vec evals_ref;
  arma::eig_sym(evals_ref, H);

  int64_t neigs = 4;

  // High-level API: lowest neigs eigenvalues with guard vectors.
  Log("lobpcg: high-level eigs_lobpcg (double)");
  EigsLobpcgResult r = eigs_lobpcg(ops, block, neigs, /*guard=*/3, /*tol=*/1e-9,
                                   /*max_iterations=*/500);
  REQUIRE(r.eigenvalues.size() == neigs);
  for (int64_t i = 0; i < neigs; ++i) {
    REQUIRE(std::abs(r.eigenvalues(i) - evals_ref(i)) < 1e-6);
  }
  // Residuals of the returned eigenpairs must be small.
  REQUIRE(r.residual_norms.max() < 1e-5);
  // History has one column per block vector and at least one row.
  REQUIRE(r.eigenvalue_history.n_cols == neigs + 3);
  REQUIRE(r.eigenvalue_history.n_rows >= 1);

  // Low-level core with matrix-block operations.
  Log("lobpcg: low-level core (double)");
  auto mult = [&](arma::mat const &X, arma::mat &Y) {
    apply(ops, block, X, block, Y);
  };
  auto dot = [&](arma::mat const &V, arma::mat const &W) {
    return math::matrix_dot(block, V, W);
  };
  int64_t nrows = size(block);
  int64_t blocksize = neigs + 3;
  arma::mat X0(nrows, blocksize, arma::fill::randn);
  linalg::lobpcg_result_t rc = linalg::lobpcg(mult, dot, X0, neigs, 1e-9, 500);
  REQUIRE(rc.criterion == "converged");
  for (int64_t i = 0; i < neigs; ++i) {
    REQUIRE(std::abs(rc.eigenvalues(i) - evals_ref(i)) < 1e-6);
  }
}

TEST_CASE("lobpcg_degeneracy", "[linalg]") {
  // An SU(2)-symmetric Heisenberg chain in the FULL spin Hilbert space (all Sz
  // sectors) has energy levels that are total-spin multiplets, degenerate with
  // multiplicity 2S+1 = 1, 3, 5, 7, ... . (In a fixed-Sz sector each multiplet
  // contributes only one state, so these degeneracies would be invisible.)
  // Unlike single-vector Lanczos, the LOBPCG block method should return the
  // lowest levels with the correct multiplicity, provided the block size
  // (neigs + guard) is large enough to hold a multiplet straddling the
  // neigs-th eigenvalue.
  int nsites = 8;
  OpSum ops;
  ops["J"] = 1.0;
  for (int s = 0; s + 1 < nsites; ++s) { // open boundary conditions
    ops += "J" * Op("SdotS", {s, s + 1});
  }
  auto block = Spinhalf(nsites); // full 2^nsites Hilbert space
  int64_t d = dim(block);

  // Dense reference spectrum and its degeneracy pattern.
  arma::mat H = matrix(ops, block);
  arma::vec evals_ref;
  arma::eig_sym(evals_ref, H);
  std::vector<int64_t> mult_ref = multiplicities(evals_ref, 1e-8);

  Log("Heisenberg N={} (full space): lowest multiplicities {} {} {} {} ...",
      nsites, mult_ref[0], mult_ref[1], mult_ref[2], mult_ref[3]);

  // Cover the first three distinct energy levels; the guard covers the fourth
  // multiplet so a degenerate cluster at the boundary is protected.
  int64_t ncover = 3;
  int64_t neigs = 0;
  for (int64_t k = 0; k < ncover; ++k) {
    REQUIRE(mult_ref[k] % 2 == 1); // SU(2): every level is (2S+1)-fold
    neigs += mult_ref[k];
  }
  int64_t guard = mult_ref[ncover] + 2;
  REQUIRE(5 * (neigs + guard) <= d); // LOBPCG needs a large enough n/blocksize

  EigsLobpcgResult r = eigs_lobpcg(ops, block, neigs, guard, 1e-9, 2000);

  
  // The lowest neigs eigenvalues (including every degenerate copy) match dense.
  for (int64_t i = 0; i < neigs; ++i) {
    REQUIRE(std::abs(r.eigenvalues(i) - evals_ref(i)) < 1e-6);
  }
  // The recovered multiplicities match the reference: no degenerate partner is
  // missed or spuriously duplicated.
  std::vector<int64_t> mult_got = multiplicities(r.eigenvalues, 1e-6);
  REQUIRE(mult_got.size() >= (size_t)ncover);
  for (int64_t k = 0; k < ncover; ++k) {
    REQUIRE(mult_got[k] == mult_ref[k]);
  }
}

TEST_CASE("lobpcg_fully_connected", "[linalg]") {
  // Uniform fully-connected Heisenberg H = sum_{i<j} S_i.S_j = 1/2(S^2 - 3N/4).
  // The energy depends only on the total spin S, so the levels form a ladder
  // E(S) = 1/2(S(S+1) - 3N/4) ordered by S, each with degeneracy D(N,S)*(2S+1)
  // where D(N,S) is the number of spin-S multiplets. The ground state (S=0 for
  // even N) is thus D(N,0)-fold degenerate -- a strong test of the block method
  // resolving a highly degenerate manifold.
  int nsites = 6;
  OpSum ops;
  ops["J"] = 1.0;
  for (int i = 0; i < nsites; ++i) {
    for (int j = i + 1; j < nsites; ++j) {
      ops += "J" * Op("SdotS", {i, j});
    }
  }
  auto block = Spinhalf(nsites); // full 2^nsites Hilbert space
  int64_t d = dim(block);

  arma::mat H = matrix(ops, block);
  arma::vec evals_ref;
  arma::eig_sym(evals_ref, H);
  std::vector<int64_t> mult_ref = multiplicities(evals_ref, 1e-8);
  Log("Fully-connected Heisenberg N={}: E(S) ladder degeneracies {} {} {} {} "
      "(= D(N,S)*(2S+1))",
      nsites, mult_ref[0], mult_ref[1], mult_ref[2], mult_ref[3]);

  // Analytic ground energy E(S=0) = -3N/8 and its D(N,0)-fold degeneracy.
  double e0_analytic = -3.0 * nsites / 8.0;
  REQUIRE(std::abs(evals_ref(0) - e0_analytic) < 1e-10);

  int64_t neigs = mult_ref[0]; // the full degenerate ground manifold
  int64_t guard = 3;
  REQUIRE(5 * (neigs + guard) <= d);

  EigsLobpcgResult r = eigs_lobpcg(ops, block, neigs, guard, 1e-9, 2000);
  for (int64_t i = 0; i < neigs; ++i) {
    REQUIRE(std::abs(r.eigenvalues(i) - evals_ref(i)) < 1e-6);
    REQUIRE(std::abs(r.eigenvalues(i) - e0_analytic) < 1e-6);
  }
  std::vector<int64_t> mult_got = multiplicities(r.eigenvalues, 1e-6);
  REQUIRE(mult_got.size() == 1);
  REQUIRE(mult_got[0] == neigs); // whole ground manifold recovered
}

TEST_CASE("lobpcg_coverage", "[linalg]") {
  using namespace xdiag::testcases::electron;
  int nsites = 6;
  OpSum ops = freefermion_alltoall(nsites);
  auto block = Electron(nsites, 3, 3);

  arma::mat H = matrix(ops, block);
  arma::vec evals_ref;
  arma::eig_sym(evals_ref, H);

  // Single eigenvalue with no guard vectors (block size == neigs == 1).
  {
    EigsLobpcgResult r = eigs_lobpcg(ops, block, 1, /*guard=*/0, 1e-10, 800);
    REQUIRE(r.criterion == "converged");
    REQUIRE(std::abs(r.eigenvalues(0) - evals_ref(0)) < 1e-8);
  }

  // Eigenvector correctness: returned vectors satisfy A v = lambda v and are
  // orthonormal.
  {
    int64_t neigs = 4;
    EigsLobpcgResult r = eigs_lobpcg(ops, block, neigs, 3, 1e-10, 1000);
    arma::mat V = r.eigenvectors.matrix();
    arma::mat AV(V.n_rows, V.n_cols);
    apply(ops, block, V, block, AV);
    for (int64_t i = 0; i < neigs; ++i) {
      REQUIRE(arma::norm(AV.col(i) - r.eigenvalues(i) * V.col(i)) < 1e-6);
    }
    arma::mat gram = V.t() * V;
    REQUIRE(arma::norm(gram - arma::eye(neigs, neigs), "inf") < 1e-8);
  }

  // Forced non-convergence: an unreachable tolerance in few iterations returns
  // criterion "maxiterations" with the best (finite) iterate.
  {
    EigsLobpcgResult r = eigs_lobpcg(ops, block, 3, 3, 1e-14, /*maxiter=*/2);
    REQUIRE(r.criterion == "maxiterations");
    REQUIRE(r.eigenvalues.size() == 3);
    REQUIRE(std::isfinite(r.residual_norms.max()));
  }
}

TEST_CASE("lobpcg_complex", "[linalg]") {
  // Complex-Hermitian coverage: a translation-symmetric Heisenberg chain at a
  // non-real momentum sector (k != 0, N/2) has a complex Hermitian matrix, so
  // eigs_lobpcg takes the complex code path (matrixC, complex eig_sym_gen).
  using namespace xdiag::testcases::spinhalf;
  int nsites = 12;
  int nup = 6;
  int k = 1; // complex momentum
  OpSum ops = HBchain(nsites, 1.0);
  auto irrep = cyclic_group_irrep(nsites, k);
  auto block = Spinhalf(nsites, nup, irrep);
  REQUIRE_FALSE(isreal(block)); // ensure we exercise the complex path
  int64_t d = dim(block);

  arma::cx_mat H = matrixC(ops, block);
  arma::vec evals_ref;
  arma::eig_sym(evals_ref, H);

  int64_t neigs = 2;
  int64_t guard = 3;
  REQUIRE(5 * (neigs + guard) <= d);

  EigsLobpcgResult r = eigs_lobpcg(ops, block, neigs, guard, 1e-9, 1000);
  for (int64_t i = 0; i < neigs; ++i) {
    REQUIRE(std::abs(r.eigenvalues(i) - evals_ref(i)) < 1e-6);
  }
  // Complex eigenvectors satisfy A v = lambda v.
  arma::cx_mat V = r.eigenvectors.matrixC();
  arma::cx_mat AV(V.n_rows, V.n_cols);
  apply(ops, block, V, block, AV);
  for (int64_t i = 0; i < neigs; ++i) {
    REQUIRE(arma::norm(AV.col(i) - r.eigenvalues(i) * V.col(i)) < 1e-6);
  }
}

TEST_CASE("lobpcg_csr", "[linalg]") {
  // Sparse-matrix coverage: run LOBPCG on precomputed CSR matrices, exercising
  // the CSRMatrix eigs_lobpcg overloads (real int64 and complex int64).
  using namespace xdiag::testcases::electron;
  int nsites = 6;
  OpSum ops = freefermion_alltoall(nsites);
  auto block = Electron(nsites, 3, 3);

  arma::mat Hd = matrix(ops, block);
  arma::vec evals_ref;
  arma::eig_sym(evals_ref, Hd);

  int64_t neigs = 3;
  int64_t guard = 3;

  // Real CSR (CSRMatrix<int64_t, double>).
  {
    CSRMatrix<int64_t, double> A = csr_matrix(ops, block);
    A.ishermitian = true;
    EigsLobpcgResult r = eigs_lobpcg(A, block, neigs, guard, 1e-9, 800);
    for (int64_t i = 0; i < neigs; ++i) {
      REQUIRE(std::abs(r.eigenvalues(i) - evals_ref(i)) < 1e-6);
    }
  }
  // Complex CSR (CSRMatrix<int64_t, complex>) -> complex code path.
  {
    CSRMatrix<int64_t, complex> A = csr_matrixC(ops, block);
    A.ishermitian = true;
    EigsLobpcgResult r = eigs_lobpcg(A, block, neigs, guard, 1e-9, 800);
    for (int64_t i = 0; i < neigs; ++i) {
      REQUIRE(std::abs(r.eigenvalues(i) - evals_ref(i)) < 1e-6);
    }
  }
}
