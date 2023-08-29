#include "eigs_lanczos.h"

#include <hydra/algebra/apply.h>
#include <hydra/algorithms/lanczos/eigvals_lanczos.h>
#include <hydra/algorithms/lanczos/lanczos.h>
#include <hydra/algorithms/lanczos/lanczos_convergence.h>

#include <hydra/states/random_state.h>
#include <hydra/utils/timing.h>

namespace hydra {

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
    HydraThrow(std::runtime_error, "Error diagonalizing tridiagonal matrix");
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
    auto dot = [](arma::cx_vec const &v, arma::cx_vec const &w) {
      return arma::cdot(v, w);
    };
    auto operation = [&eigenvectors, &revecs, &iter,
                      neigvals](arma::cx_vec const &v) {
      eigenvectors.matrixC(false) +=
          kron(v, revecs.submat(iter - 1, 0, iter - 1, neigvals - 1));
    };

    lanczos::lanczos(mult, dot, converged, operation, v0, r.niterations,
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
    auto dot = [](arma::vec const &v, arma::vec const &w) {
      return arma::dot(v, w);
    };
    auto operation = [&eigenvectors, &revecs, &iter,
                      neigvals](arma::vec const &v) {
      eigenvectors.matrix(false) +=
          kron(v, revecs.submat(iter - 1, 0, iter - 1, neigvals - 1));
    };
    lanczos::lanczos(mult, dot, converged, operation, v0, r.niterations,
                     deflation_tol);
  }

  return {r.alphas,     r.betas,       r.eigenvalues,
          eigenvectors, r.niterations, r.criterion};
} catch (...) {
  HydraRethrow("Error performing eigenvector Lanczos algorithm");
  return eigs_lanczos_result_t();
}
} // namespace hydra
