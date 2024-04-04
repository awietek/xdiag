#include "eigs_lanczos.h"
#include <xdiag/algebra/algebra.h>
#include <xdiag/algebra/apply.h>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.h>
#include <xdiag/algorithms/lanczos/lanczos.h>
#include <xdiag/algorithms/lanczos/lanczos_convergence.h>

#include <xdiag/states/random_state.h>
#include <xdiag/utils/print_macro.h>
#include <xdiag/utils/timing.h>

namespace xdiag {

eigs_lanczos_result_t eigs_lanczos(BondList const &bonds,
                                   block_variant_t const &block,
                                   int64_t neigvals, double precision,
                                   int64_t max_iterations, bool force_complex,
                                   double deflation_tol,
                                   int64_t random_seed) try {
  if (neigvals < 1) {
    throw(std::invalid_argument("Argument \"neigvals\" needs to be >= 1"));
  }
  if (!bonds.ishermitian()) {
    throw(std::invalid_argument("Input BondList is not hermitian"));
  }

  bool cplx = bonds.iscomplex() || iscomplex(block) || force_complex;
  State state0(block, !cplx);
  fill(state0, RandomState(random_seed));

  // Perform first run to compute eigenvalues
  auto r = eigvals_lanczos(bonds, block, neigvals, precision, max_iterations,
                           force_complex, deflation_tol, random_seed);

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
    XDiagThrow(std::runtime_error, "Error diagonalizing tridiagonal matrix");
  }

  auto converged = [](Tmatrix const &) -> bool { return false; };

  fill(state0, RandomState(random_seed));
  State eigenvectors(block, !cplx, neigvals);

  int64_t iter = 1;
  // Setup complex Lanczos run
  if (cplx) {
    arma::cx_vec v0 = state0.vectorC(0, false);
    auto mult = [&iter, &bonds, &block](arma::cx_vec const &v,
                                        arma::cx_vec &w) {
      auto ta = rightnow();
      apply(bonds, block, v, block, w);
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
    arma::vec v0 = state0.vector(0, false);
    auto mult = [&iter, &bonds, &block](arma::vec const &v, arma::vec &w) {
      auto ta = rightnow();
      apply(bonds, block, v, block, w);
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
} catch (...) {
  XDiagRethrow("Error performing eigenvector Lanczos algorithm");
  return eigs_lanczos_result_t();
}
} // namespace xdiag
