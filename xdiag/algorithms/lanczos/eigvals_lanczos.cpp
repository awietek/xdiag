#include "eigvals_lanczos.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algorithms/lanczos/lanczos.hpp>
#include <xdiag/algorithms/lanczos/lanczos_convergence.hpp>

#include <xdiag/states/random_state.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/utils/timing.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/parallel/mpi/cdot_distributed.hpp>
#endif

namespace xdiag {

eigvals_lanczos_result_t
eigvals_lanczos(OpSum const &ops, Block const &block, State &state0,
                int64_t neigvals, double precision, int64_t max_iterations,
                bool force_complex, double deflation_tol) try {
  if (neigvals < 1) {
    XDIAG_THROW("Argument \"neigvals\" needs to be >= 1");
  }
  // if (!ops.ishermitian()) {
  //   XDIAG_THROW("Input OpSum is not hermitian");
  // }
  bool cplx =
      !ops.isreal() || !isreal(block) || force_complex || !state0.isreal();
  if (!state0.isreal() && !isreal(block)) {
    Log(1,
        "warning: starting REAL block diagonalization with COMPLEX startstate");
  }
  auto converged = [neigvals, precision](Tmatrix const &tmat) -> bool {
    return lanczos::converged_eigenvalues(tmat, neigvals, precision);
  };
  lanczos::lanczos_result_t r;
  int64_t iter = 1;
  // Setup complex Lanczos run
  if (cplx) {
    state0.make_complex();
    arma::cx_vec v0 = state0.vectorC(0, false);
    auto mult = [&iter, &ops, &block](arma::cx_vec const &v, arma::cx_vec &w) {
      auto ta = rightnow();
      apply(ops, block, v, block, w);
      Log(1, "Lanczos iteration {}", iter);
      timing(ta, rightnow(), "MVM", 1);
      ++iter;
    };

    auto operation = [](arma::cx_vec const &) {};
    auto dotf = [&block](arma::cx_vec const &v, arma::cx_vec const &w) {
      return dot(block, v, w);
    };
    r = lanczos::lanczos(mult, dotf, converged, operation, v0, max_iterations,
                         deflation_tol);

    // Setup real Lanczos run
  } else {
    arma::vec v0 = state0.vector(0, false);
    auto mult = [&iter, &ops, &block](arma::vec const &v, arma::vec &w) {
      auto ta = rightnow();
      apply(ops, block, v, block, w);
      Log(1, "Lanczos iteration {}", iter);
      timing(ta, rightnow(), "MVM", 1);
      ++iter;
    };

    auto operation = [](arma::vec const &) {};
    auto dotf = [&block](arma::vec const &v, arma::vec const &w) {
      return dot(block, v, w);
    };
    r = lanczos::lanczos(mult, dotf, converged, operation, v0, max_iterations,
                         deflation_tol);
  }
  return {r.alphas, r.betas, r.eigenvalues, r.niterations, r.criterion};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return eigvals_lanczos_result_t();
}

eigvals_lanczos_result_t
eigvals_lanczos(OpSum const &ops, Block const &block, int64_t neigvals,
                double precision, int64_t max_iterations, bool force_complex,
                double deflation_tol, int64_t random_seed) try {

  if (neigvals < 1) {
    XDIAG_THROW("Argument \"neigvals\" needs to be >= 1");
  }
  // if (!ops.ishermitian()) {
  //   XDIAG_THROW("Input OpSum is not hermitian");
  // }

  bool cplx = (!ops.isreal()) || !isreal(block) || force_complex;
  State state0(block, !cplx);
  fill(state0, RandomState(random_seed));

  auto r = eigvals_lanczos(ops, block, state0, neigvals, precision,
                           max_iterations, force_complex, deflation_tol);

  return {r.alphas, r.betas, r.eigenvalues, r.niterations, r.criterion};

} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return eigvals_lanczos_result_t();
}
} // namespace xdiag
