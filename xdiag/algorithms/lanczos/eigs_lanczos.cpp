#include "eigs_lanczos.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
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

EigsLanczosResult eigs_lanczos(OpSum const &ops, State const &state0,
                               int64_t neigvals, double precision,
                               int64_t max_iterations,
                               double deflation_tol) try {
  if (neigvals < 1) {
    XDIAG_THROW("Argument \"neigvals\" needs to be >= 1");
  }
  if (!isapprox(ops, hc(ops))) {
    XDIAG_THROW("Input OpSum is not hermitian");
  }
  auto const &block = state0.block();

  bool real = isreal(ops) && isreal(block) && isreal(state0);

  State state1 = state0;
  if (!real) {
    state1.make_complex();
  }
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

  auto converged = [](Tmatrix const &) -> bool { return false; };

  State eigenvectors(block, real, neigvals);
  state1 = state0;
  if (!real) {
    state1.make_complex();
  }

  int64_t iter = 1;
  // Setup complex Lanczos run
  if (!real) {
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

    // Setup real Lanczos run
  } else {
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
  }

  return {r.alphas,     r.betas,       r.eigenvalues,
          eigenvectors, r.niterations, r.criterion};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

// starting from random vector
EigsLanczosResult eigs_lanczos(OpSum const &ops, Block const &block,
                               int64_t neigvals, double precision,
                               int64_t max_iterations, double deflation_tol,
                               int64_t random_seed) try {
  if (neigvals < 1) {
    XDIAG_THROW("Argument \"neigvals\" needs to be >= 1");
  }
  // if (!ops.ishermitian()) {
  //   XDIAG_THROW("Input OpSum is not hermitian");
  // }

  bool real = isreal(ops) && isreal(block);
  State state0(block, real);
  fill(state0, RandomState(random_seed));

  auto r = eigs_lanczos(ops, state0, neigvals, precision, max_iterations,
                        deflation_tol);

  return {r.alphas,       r.betas,       r.eigenvalues,
          r.eigenvectors, r.niterations, r.criterion};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
