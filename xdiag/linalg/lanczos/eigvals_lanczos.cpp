// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "eigvals_lanczos.hpp"

#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/sparse/apply.hpp>
#include <xdiag/linalg/lanczos/lanczos.hpp>
#include <xdiag/linalg/lanczos/lanczos_convergence.hpp>
#include <xdiag/math/dot.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/timing.hpp>

#ifdef XDIAG_DISTRIBUTED
#include <xdiag/mpi/cdot_distributed.hpp>
#endif

namespace xdiag {

// Sets up the multiply/dot/operation callbacks for a given coefficient type
// and runs the bare Lanczos iteration. Shared by the real (coeff_t = double)
// and complex (coeff_t = complex) code paths.
template <typename coeff_t, typename op_t, typename converged_f>
static lanczos::lanczos_result_t
run_eigvals_lanczos(op_t const &ops, Block const &block, arma::Col<coeff_t> &v0,
                    converged_f converged, int64_t max_iterations,
                    double deflation_tol) {
  int64_t iter = 1;
  auto mult = [&](arma::Col<coeff_t> const &v, arma::Col<coeff_t> &w) {
    auto ta = rightnow();
    apply(ops, block, v, block, w);
    Log(1, "Lanczos iteration {}", iter);
    timing(ta, rightnow(), "MVM", 1);
    ++iter;
  };
  auto operation = [](arma::Col<coeff_t> const &) {};
  // different dot product for normal and distributed blocks
  auto dotf = [&](arma::Col<coeff_t> const &v, arma::Col<coeff_t> const &w) {
    return math::dot(block, v, w);
  };
  return lanczos::lanczos(mult, dotf, converged, operation, v0, max_iterations,
                          deflation_tol);
}

///////////////////////////////////////////////////////////////
// Routine with random state initialization

template <typename op_t>
static EigvalsLanczosResult
eigvals_lanczos(op_t const &ops, Block const &block, int64_t neigvals,
                double precision, int64_t max_iterations, double deflation_tol,
                int64_t random_seed) try {
  bool real = isreal(ops) && isreal(block);
  State state0(block, real);
  fill(state0, RandomState(random_seed));
  return eigvals_lanczos_inplace(ops, state0, neigvals, precision,
                                 max_iterations, deflation_tol);
}
XDIAG_CATCH

EigvalsLanczosResult eigvals_lanczos(OpSum const &ops, Block const &block,
                                     int64_t neigvals, double precision,
                                     int64_t max_iterations,
                                     double deflation_tol,
                                     int64_t random_seed) try {
  return eigvals_lanczos<OpSum>(ops, block, neigvals, precision, max_iterations,
                                deflation_tol, random_seed);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
EigvalsLanczosResult
eigvals_lanczos(CSRMatrix<idx_t, coeff_t> const &ops, Block const &block,
                int64_t neigvals, double precision, int64_t max_iterations,
                double deflation_tol, int64_t random_seed) try {
  return eigvals_lanczos<CSRMatrix<idx_t, coeff_t>>(ops, block, neigvals,
                                                    precision, max_iterations,
                                                    deflation_tol, random_seed);
}
XDIAG_CATCH

///////////////////////////////////////////////////////////////
// Routine with given starting state which is copied

EigvalsLanczosResult eigvals_lanczos(OpSum const &ops, State psi0,
                                     int64_t neigvals, double precision,
                                     int64_t max_iterations,
                                     double deflation_tol) try {
  return eigvals_lanczos_inplace(ops, psi0, neigvals, precision, max_iterations,
                                 deflation_tol);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
EigvalsLanczosResult eigvals_lanczos(CSRMatrix<idx_t, coeff_t> const &ops,
                                     State psi0, int64_t neigvals,
                                     double precision, int64_t max_iterations,
                                     double deflation_tol) try {
  return eigvals_lanczos_inplace(ops, psi0, neigvals, precision, max_iterations,
                                 deflation_tol);
}
XDIAG_CATCH

///////////////////////////////////////////////////////////////
// Routine with given starting state which is overwritten

template <typename op_t>
static EigvalsLanczosResult
eigvals_lanczos_inplace(op_t const &ops, State &psi0, int64_t neigvals,
                        double precision, int64_t max_iterations,
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

  if (!ishermitian(ops, psi0.block())) {
    XDIAG_THROW("Input operator is not Hermitian. The Lanczos algorithm can "
                "only be applied to Hermitian operators.");
  }

  auto const &block = psi0.block();
  bool real = isreal(ops) && isreal(block) && isreal(psi0);

  auto converged = [neigvals, precision](Tmatrix const &tmat) -> bool {
    return lanczos::converged_eigenvalues(tmat, neigvals, precision);
  };

  lanczos::lanczos_result_t r;
  if (real) {                             // Real Lanczos algorithm
    arma::vec v0 = psi0.vector(0, false); // not copied
    r = run_eigvals_lanczos(ops, block, v0, converged, max_iterations,
                            deflation_tol);
  } else { // Complex Lanczos algorithm
    psi0.make_complex();
    arma::cx_vec v0 = psi0.vectorC(0, false); // not copied
    r = run_eigvals_lanczos(ops, block, v0, converged, max_iterations,
                            deflation_tol);
  }
  return {r.alphas, r.betas, r.eigenvalues, r.niterations, r.criterion};
}
XDIAG_CATCH

EigvalsLanczosResult eigvals_lanczos_inplace(OpSum const &ops, State &psi0,
                                             int64_t neigvals, double precision,
                                             int64_t max_iterations,
                                             double deflation_tol) try {
  return eigvals_lanczos_inplace<OpSum>(ops, psi0, neigvals, precision,
                                        max_iterations, deflation_tol);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
EigvalsLanczosResult
eigvals_lanczos_inplace(CSRMatrix<idx_t, coeff_t> const &ops, State &psi0,
                        int64_t neigvals, double precision,
                        int64_t max_iterations, double deflation_tol) try {
  return eigvals_lanczos_inplace<CSRMatrix<idx_t, coeff_t>>(
      ops, psi0, neigvals, precision, max_iterations, deflation_tol);
}
XDIAG_CATCH

// Template instantiations for every (idx_t, coeff_t) sparse-matrix combination
#define XDIAG_INST(IDX, COEFF)                                                 \
  template EigvalsLanczosResult eigvals_lanczos(                               \
      CSRMatrix<IDX, COEFF> const &, Block const &, int64_t, double, int64_t,  \
      double, int64_t);                                                        \
  template EigvalsLanczosResult eigvals_lanczos(                               \
      CSRMatrix<IDX, COEFF> const &, State, int64_t, double, int64_t, double); \
  template EigvalsLanczosResult eigvals_lanczos_inplace(                       \
      CSRMatrix<IDX, COEFF> const &, State &, int64_t, double, int64_t,        \
      double);

XDIAG_INST(int32_t, double)
XDIAG_INST(int32_t, complex)
XDIAG_INST(int64_t, double)
XDIAG_INST(int64_t, complex)
#undef XDIAG_INST

} // namespace xdiag
