// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "eigvals_lanczos.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algorithms/lanczos/lanczos.hpp>
#include <xdiag/algorithms/lanczos/lanczos_convergence.hpp>

#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>

#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/real.hpp>

#include <xdiag/utils/timing.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/parallel/mpi/cdot_distributed.hpp>
#endif

namespace xdiag {

EigvalsLanczosResult eigvals_lanczos(OpSum const &ops, Block const &block,
                                     int64_t neigvals, double precision,
                                     int64_t max_iterations,
                                     double deflation_tol,
                                     int64_t random_seed) try {

  if (neigvals < 1) {
    XDIAG_THROW("Argument \"neigvals\" needs to be >= 1");
  } else if (neigvals > dim(block)) {
    neigvals = dim(block);
  }

  bool real = isreal(ops) && isreal(block);
  State state0(block, real);
  fill(state0, RandomState(random_seed));

  auto r = eigvals_lanczos_inplace(ops, state0, neigvals, precision,
                                   max_iterations, deflation_tol);

  return {r.alphas, r.betas, r.eigenvalues, r.niterations, r.criterion};

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

EigvalsLanczosResult eigvals_lanczos(OpSum const &ops, State psi0,
                                     int64_t neigvals, double precision,
                                     int64_t max_iterations,
                                     double deflation_tol) try {
  return eigvals_lanczos_inplace(ops, psi0, neigvals, precision, max_iterations,
                                 deflation_tol);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

EigvalsLanczosResult eigvals_lanczos_inplace(OpSum const &ops, State &psi0,
                                             int64_t neigvals, double precision,
                                             int64_t max_iterations,
                                             double deflation_tol) try {
  if (dim(psi0) == 0) {
    Log.warn(
        "Warning: initial state zero dimensional in eigvals_lanczos_inplace");
    return EigvalsLanczosResult();
  }

  if (neigvals < 1) {
    XDIAG_THROW("Argument \"neigvals\" needs to be >= 1");
  } else if (neigvals > dim(psi0.block())) {
    neigvals = dim(psi0.block());
  }

  if (!isvalid(psi0)) {
    XDIAG_THROW("Initial state must be a valid state (i.e. not default "
                "constructed by e.g. an annihilation operator)");
  }

  if (!isapprox(ops, hc(ops))) {
    XDIAG_THROW("Input OpSum is not hermitian");
  }

  auto const &block = psi0.block();

  bool real = isreal(ops) && isreal(block) && isreal(psi0);
  auto converged = [neigvals, precision](Tmatrix const &tmat) -> bool {
    return lanczos::converged_eigenvalues(tmat, neigvals, precision);
  };
  lanczos::lanczos_result_t r;
  int64_t iter = 1;
  // Setup complex Lanczos run
  if (!real) {
    psi0.make_complex();
    arma::cx_vec v0 = psi0.vectorC(0, false);
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
    arma::vec v0 = psi0.vector(0, false);
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
}

} // namespace xdiag
