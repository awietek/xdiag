#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/common.h>
#include <hydra/linalg/lanczos/tmatrix.h>
#include <hydra/utils/logger.h>

#include <hydra/linalg/lanczos/lanczos_step.h>

namespace hydra {

// Generic Lanczos implementation building multiple vectors
template <class coeff_t, class multiply_f, class convergence_f>
Tmatrix lanczos(multiply_f mult, arma::Col<coeff_t> &v0,
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

    lanczos_step(v0, v1, w, alpha, beta, mult, dot);
    tmatrix.append(alpha, beta);
    tmatrix.print_log();

    // Finish if Lanczos sequence is exhausted
    if (std::abs(beta) > deflation_tol) {
      v1 /= beta;
    } else {
      break;
    }

    ++iteration;
    if (iteration >= max_iterations)
      break;
  }

  return tmatrix;
}

} // namespace hydra
