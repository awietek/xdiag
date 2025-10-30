// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "eigs_lanczos.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/sparse/apply.hpp>
#include <xdiag/algebra/sparse/logic.hpp>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/algorithms/lanczos/lanczos.hpp>
#include <xdiag/algorithms/lanczos/lanczos_convergence.hpp>

#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>

#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/real.hpp>

#include <xdiag/utils/timing.hpp>

namespace xdiag {


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
  if (!isapprox(ops, hc(ops))) {
    XDIAG_THROW("Input OpSum is not hermitian");
  }
  auto const &block = state0.block();

  bool real = isreal(ops) && isreal(block) && isreal(state0);


  // store initial state, such that it can be used again in second run
  State state1 = state0;

  // Perform first run to compute eigenvalues
  auto r = eigvals_lanczos_inplace(ops, state1, neigvals, precision,
                                   max_iterations, deflation_tol);

  // Perform second run to compute the eigenvectors
  arma::mat tmat = arma::diagmat(r.alphas);
  if (r.alphas.n_rows > 1) {
    tmat += arma::diagmat(r.betas.head(r.betas.size() - 1), 1) +
            arma::diagmat(r.betas.head(r.betas.size() - 1), -1);
  }

  arma::vec reigs;
  arma::mat revecs;
  try {
    arma::eig_sym(reigs, revecs, tmat);
  } catch (...) {
    XDIAG_THROW("Error diagonalizing tridiagonal matrix");
  }

  // now we prepare second run, no convergence is checked, just fixed number
  // of iterations is performed
  auto converged = [](Tmatrix const &) -> bool { return false; };
  auto const &block = state0.block();
  State eigenvectors(block, isreal(ops), neigvals);
  state1 = state0;
  int64_t iter = 1;

  if (isreal(ops) && isreal(block) && isreal(state1)) { // Real Lanczos
    arma::vec v0 = state1.vector(0, false);
    auto mult = [&iter, &ops, &block](arma::vec const &v, arma::vec &w) {
      auto ta = rightnow();
      apply(ops, block, v, block, w);
      Log(1, "Lanczos iteration {}", iter);
      timing(ta, rightnow(), "MVM", 1);
      ++iter;
    };
    auto dotf = [&block](arma::vec const &v, arma::vec const &w) {
      return dot(block, v, w);
    };
    auto operation = [&eigenvectors, &revecs, &iter,
                      neigvals](arma::vec const &v) {
      eigenvectors.matrix(false) +=
          kron(v, revecs.submat(iter - 1, 0, iter - 1, neigvals - 1));
    };
    lanczos::lanczos(mult, dotf, converged, operation, v0, r.niterations,
                     deflation_tol);

  } else { // Complex Lanczos
    arma::cx_vec v0 = state1.vectorC(0, false);
    auto mult = [&iter, &ops, &block](arma::cx_vec const &v, arma::cx_vec &w) {
      auto ta = rightnow();
      apply(ops, block, v, block, w);
      Log(1, "Lanczos iteration (rerun) {}", iter);
      timing(ta, rightnow(), "MVM", 1);
      ++iter;
    };
    auto dotf = [&block](arma::cx_vec const &v, arma::cx_vec const &w) {
      return dot(block, v, w);
    };
    auto operation = [&eigenvectors, &revecs, &iter,
                      neigvals](arma::cx_vec const &v) {
      eigenvectors.matrixC(false) +=
          kron(v, revecs.submat(iter - 1, 0, iter - 1, neigvals - 1));
    };

    lanczos::lanczos(mult, dotf, converged, operation, v0, r.niterations,
                     deflation_tol);
  }

  return {r.alphas,     r.betas,       r.eigenvalues,
          eigenvectors, r.niterations, r.criterion};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

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
} XDIAG_CATCH

EigsLanczosResult eigs_lanczos(OpSum const &ops, Block const &block,
                               int64_t neigvals, double precision,
                               int64_t max_iterations, double deflation_tol,
                               int64_t random_seed) try {
  return eigs_lanczos<OpSum>(ops, block, neigvals, precision, max_iterations,
                             deflation_tol, random_seed);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
EigsLanczosResult eigs_lanczos(CSRMatrix<idx_t, coeff_t> const &ops,
                               Block const &block, int64_t neigvals,
                               double precision, int64_t max_iterations,
                               double deflation_tol, int64_t random_seed) try {
  return eigs_lanczos<CSRMatrix<idx_t, coeff_t>>(ops, block, neigvals,
                                                 precision, max_iterations,
                                                 deflation_tol, random_seed);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template EigsLanczosResult
eigs_lanczos(CSRMatrix<int32_t, double> const &ops, Block const &block,
             int64_t neigvals, double precision, int64_t max_iterations,
             double deflation_tol, int64_t random_seed);
template EigsLanczosResult
eigs_lanczos(CSRMatrix<int32_t, complex> const &ops, Block const &block,
             int64_t neigvals, double precision, int64_t max_iterations,
             double deflation_tol, int64_t random_seed);
template EigsLanczosResult
eigs_lanczos(CSRMatrix<int64_t, double> const &ops, Block const &block,
             int64_t neigvals, double precision, int64_t max_iterations,
             double deflation_tol, int64_t random_seed);
template EigsLanczosResult
eigs_lanczos(CSRMatrix<int64_t, complex> const &ops, Block const &block,
             int64_t neigvals, double precision, int64_t max_iterations,
             double deflation_tol, int64_t random_seed);

///////////////////////////////////////////////////////////////
// Routine with random state initialization
EigsLanczosResult eigs_lanczos(OpSum const &ops, State const &state0,
                               int64_t neigvals, double precision,
                               int64_t max_iterations,
                               double deflation_tol) try {
  return eigs_lanczos<OpSum>(ops, state0, neigvals, precision, max_iterations,
                             deflation_tol);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
EigsLanczosResult eigs_lanczos(CSRMatrix<idx_t, coeff_t> const &ops,
                               State const &state0, int64_t neigvals,
                               double precision, int64_t max_iterations,
                               double deflation_tol) try {
  return eigs_lanczos<CSRMatrix<idx_t, coeff_t>>(
      ops, state0, neigvals, precision, max_iterations, deflation_tol);
} XDIAG_CATCH

template EigsLanczosResult eigs_lanczos(CSRMatrix<int32_t, double> const &ops,
                                        State const &state0, int64_t neigvals,
                                        double precision,
                                        int64_t max_iterations,
                                        double deflation_tol);
template EigsLanczosResult eigs_lanczos(CSRMatrix<int32_t, complex> const &ops,
                                        State const &state0, int64_t neigvals,
                                        double precision,
                                        int64_t max_iterations,
                                        double deflation_tol);
template EigsLanczosResult eigs_lanczos(CSRMatrix<int64_t, double> const &ops,
                                        State const &state0, int64_t neigvals,
                                        double precision,
                                        int64_t max_iterations,
                                        double deflation_tol);
template EigsLanczosResult eigs_lanczos(CSRMatrix<int64_t, complex> const &ops,
                                        State const &state0, int64_t neigvals,
                                        double precision,
                                        int64_t max_iterations,
                                        double deflation_tol);
} // namespace xdiag
