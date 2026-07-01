// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "eigs_lanczos.hpp"

#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/sparse/apply.hpp>
#include <xdiag/linalg/lanczos/eigvals_lanczos.hpp>
#include <xdiag/linalg/lanczos/lanczos.hpp>
#include <xdiag/linalg/lanczos/lanczos_convergence.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/math/dot.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag {

// Re-runs the Lanczos iteration for a given coefficient type, accumulating the
// eigenvectors on the fly via the "operation" callback. The tridiagonal
// eigenvectors "revecs" project each Lanczos vector onto the wanted Ritz
// vectors. Shared by the real (coeff_t = double) and complex code paths.
template <typename coeff_t, typename op_t>
static void run_eigs_lanczos(op_t const &ops, Block const &block,
                             arma::Col<coeff_t> &v0, arma::mat const &revecs,
                             int64_t neigvals, int64_t max_iterations,
                             double deflation_tol, State &eigenvectors) {
  int64_t iter = 1;
  auto mult = [&](arma::Col<coeff_t> const &v, arma::Col<coeff_t> &w) {
    auto ta = rightnow();
    apply(ops, block, v, block, w);
    Log(1, "Lanczos iteration {}", iter);
    timing(ta, rightnow(), "MVM", 1);
    ++iter;
  };
  auto dotf = [&](arma::Col<coeff_t> const &v, arma::Col<coeff_t> const &w) {
    return math::dot(block, v, w);
  };
  auto operation = [&](arma::Col<coeff_t> const &v) {
    auto coeffs = revecs.submat(iter - 1, 0, iter - 1, neigvals - 1);
    if constexpr (isreal<coeff_t>()) {
      eigenvectors.matrix(false) += arma::kron(v, coeffs);
    } else {
      eigenvectors.matrixC(false) += arma::kron(v, coeffs);
    }
  };
  // no convergence check: perform a fixed number of iterations
  auto converged = [](Tmatrix const &) -> bool { return false; };
  lanczos::lanczos(mult, dotf, converged, operation, v0, max_iterations,
                   deflation_tol);
}

template <typename op_t>
static EigsLanczosResult eigs_lanczos(op_t const &ops, State const &state0,
                                      int64_t neigvals, double precision,
                                      int64_t max_iterations,
                                      double deflation_tol) try {
  if (dim(state0) == 0) {
    Log.warn("Warning: initial state zero dimensional in eigs_lanczos");
    return EigsLanczosResult();
  }

  if (neigvals < 1) {
    XDIAG_THROW("Argument \"neigvals\" needs to be >= 1");
  } else if (neigvals > dim(state0.block())) {
    neigvals = dim(state0.block());
  }
  if (!isvalid(state0)) {
    XDIAG_THROW("Initial state must be a valid state (i.e. not default "
                "constructed by e.g. an annihilation operator)");
  }
  if (!ishermitian(ops, state0.block())) {
    XDIAG_THROW("Input OpSum is not hermitian");
  }

  // store initial state, such that it can be used again in second run
  State state1 = state0;

  // Perform first run to compute eigenvalues
  auto r = eigvals_lanczos_inplace(ops, state1, neigvals, precision,
                                   max_iterations, deflation_tol);

  // Perform second run to compute the eigenvectors. The tridiagonal T-matrix
  // is reconstructed from the recurrence coefficients (cf. Tmatrix::mat()).
  Tmatrix tmatrix(arma::conv_to<std::vector<double>>::from(r.alphas),
                  arma::conv_to<std::vector<double>>::from(r.betas));
  arma::mat tmat = tmatrix.mat();

  arma::vec reigs;
  arma::mat revecs;
  try {
    arma::eig_sym(reigs, revecs, tmat);
  } catch (...) {
    XDIAG_THROW("Error diagonalizing tridiagonal matrix");
  }

  // Second run: no convergence is checked, just a fixed number of iterations
  auto const &block = state0.block();
  State eigenvectors(block, isreal(ops), neigvals);
  state1 = state0;

  if (isreal(ops) && isreal(block) && isreal(state1)) { // Real Lanczos
    arma::vec v0 = state1.vector(0, false);
    run_eigs_lanczos(ops, block, v0, revecs, neigvals, r.niterations,
                     deflation_tol, eigenvectors);
  } else { // Complex Lanczos
    arma::cx_vec v0 = state1.vectorC(0, false);
    run_eigs_lanczos(ops, block, v0, revecs, neigvals, r.niterations,
                     deflation_tol, eigenvectors);
  }

  return {r.alphas,     r.betas,       r.eigenvalues,
          eigenvectors, r.niterations, r.criterion};
}
XDIAG_CATCH

template <typename op_t>
static EigsLanczosResult
eigs_lanczos(op_t const &ops, Block const &block, int64_t neigvals,
             double precision, int64_t max_iterations, double deflation_tol,
             int64_t random_seed) try {
  bool real = isreal(ops) && isreal(block);
  State state0(block, real);
  fill(state0, RandomState(random_seed));
  return eigs_lanczos(ops, state0, neigvals, precision, max_iterations,
                      deflation_tol);
}
XDIAG_CATCH

EigsLanczosResult eigs_lanczos(OpSum const &ops, Block const &block,
                               int64_t neigvals, double precision,
                               int64_t max_iterations, double deflation_tol,
                               int64_t random_seed) try {
  return eigs_lanczos<OpSum>(ops, block, neigvals, precision, max_iterations,
                             deflation_tol, random_seed);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
EigsLanczosResult eigs_lanczos(CSRMatrix<idx_t, coeff_t> const &ops,
                               Block const &block, int64_t neigvals,
                               double precision, int64_t max_iterations,
                               double deflation_tol, int64_t random_seed) try {
  return eigs_lanczos<CSRMatrix<idx_t, coeff_t>>(ops, block, neigvals,
                                                 precision, max_iterations,
                                                 deflation_tol, random_seed);
}
XDIAG_CATCH

///////////////////////////////////////////////////////////////
// Routine with random state initialization
EigsLanczosResult eigs_lanczos(OpSum const &ops, State const &state0,
                               int64_t neigvals, double precision,
                               int64_t max_iterations,
                               double deflation_tol) try {
  return eigs_lanczos<OpSum>(ops, state0, neigvals, precision, max_iterations,
                             deflation_tol);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
EigsLanczosResult eigs_lanczos(CSRMatrix<idx_t, coeff_t> const &ops,
                               State const &state0, int64_t neigvals,
                               double precision, int64_t max_iterations,
                               double deflation_tol) try {
  return eigs_lanczos<CSRMatrix<idx_t, coeff_t>>(
      ops, state0, neigvals, precision, max_iterations, deflation_tol);
}
XDIAG_CATCH

// Template instantiations for every (idx_t, coeff_t) sparse-matrix combination
#define XDIAG_INST(IDX, COEFF)                                                 \
  template EigsLanczosResult eigs_lanczos(CSRMatrix<IDX, COEFF> const &,       \
                                          Block const &, int64_t, double,      \
                                          int64_t, double, int64_t);           \
  template EigsLanczosResult eigs_lanczos(CSRMatrix<IDX, COEFF> const &,       \
                                          State const &, int64_t, double,      \
                                          int64_t, double);

XDIAG_INST(int32_t, double)
XDIAG_INST(int32_t, complex)
XDIAG_INST(int64_t, double)
XDIAG_INST(int64_t, complex)
#undef XDIAG_INST

} // namespace xdiag
