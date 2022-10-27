#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/algorithms/lanczos/lanczos_step.h>
#include <hydra/algorithms/lanczos/tmatrix.h>
#include <hydra/common.h>
#include <hydra/utils/logger.h>

namespace hydra {

// Lanczos applying a function on every Lanczos vector
template <class coeff_t, class multiply_f, class convergence_f,
          class vector_apply_f>
inline Tmatrix
lanczos_vector_apply_inplace(multiply_f mult, arma::Col<coeff_t> &v0,
                             vector_apply_f vector_apply,
                             convergence_f converged, int max_iterations = 1000,
                             double deflation_tol = 1e-7) {

#ifdef HYDRA_ENABLE_MPI
  auto dot = [](arma::Col<coeff_t> const &v,
                arma::Col<coeff_t> const &w) -> coeff_t {
    return DotMPI(v, w);
  };
#else
  auto dot = [](arma::Col<coeff_t> const &v,
                arma::Col<coeff_t> const &w) -> coeff_t {
    return arma::cdot(v, w);
  };
#endif

  auto norm = [&dot](arma::Col<coeff_t> const &v) {
    return std::sqrt(hydra::real(dot(v, v)));
  };

  auto tmatrix = Tmatrix();

  // Initialize Lanczos vectors and tmatrix
  auto v1 = v0;
  arma::Col<coeff_t> w(v0.size(), arma::fill::zeros);
  v0.zeros();
  double alpha = 0.;
  double beta = 0.;

  // Normalize start vector or return if norm is zero
  coeff_t v1_norm = norm(v1);
  if (std::abs(v1_norm) > 1e-12) {
    v1 /= v1_norm;
  } else {
    return tmatrix;
  }

  // Main Lanczos loop
  int iteration = 0;
  while (!converged(tmatrix)) {

    vector_apply(iteration, v1);

    lanczos_step(v0, v1, w, alpha, beta, mult, dot);
    tmatrix.append(alpha, beta);
    tmatrix.print_log();

    // Finish if Lanczos sequence is exhausted
    if (std::abs(beta) > deflation_tol) {
      v1 /= beta;
    } else {
      return tmatrix;
    }

    ++iteration;
    if (iteration >= max_iterations)
      break;
  }

  vector_apply(iteration, v1);

  return tmatrix;
}

template <class coeff_t, class vector_apply_f>
inline Tmatrix
lanczos_vector_apply_inplace(BondList const &bonds, State<coeff_t> &state_0,
                             vector_apply_f vector_apply, int max_iterations,
                             double deflation_tol = 1e-7) {
  auto const &block = state_0.block();
  auto &v0 = state_0.vector();

  int iter = 1;
  auto mult = [&iter, &bonds, &block](arma::Col<coeff_t> const &v,
                                      arma::Col<coeff_t> &w) {
    auto ta = rightnow();
    apply(bonds, block, v, block, w);
    Log(1, "Lanczos iteration {}", iter);
    timing(ta, rightnow(), "MVM", 1);
    ++iter;
  };

  auto t0 = rightnow();
  auto tmat = lanczos_vector_apply_inplace(
      mult, v0, vector_apply, [](Tmatrix const &) { return false; },
      max_iterations, deflation_tol);
  timing(t0, rightnow(), "Lanczos time", 1);

  return tmat;
}

template <class coeff_t, class vector_apply_f>
Tmatrix lanczos_vector_apply(BondList const &bonds, State<coeff_t> state_0,
                             vector_apply_f vector_apply, int max_iterations,
                             double deflation_tol = 1e-7) {
  return lanczos_vector_apply_inplace(bonds, state_0, vector_apply,
                                      max_iterations, deflation_tol);
}

} // namespace hydra
