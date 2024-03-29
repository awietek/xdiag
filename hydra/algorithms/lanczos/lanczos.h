#pragma once

#include <hydra/extern/armadillo/armadillo>

#include <hydra/algorithms/lanczos/lanczos_step.h>
#include <hydra/algorithms/lanczos/tmatrix.h>
#include <hydra/common.h>
#include <hydra/utils/logger.h>

namespace hydra::lanczos {

struct lanczos_result_t {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  int64_t niterations;
  std::string criterion;
};

template <class coeff_t, class mult_f, class dot_f, class converged_f,
          class operation_f>
lanczos_result_t lanczos(mult_f mult, dot_f dot, converged_f converged,
                         operation_f operation, arma::Col<coeff_t> &v0,
                         int max_iterations = 1000,
                         double deflation_tol = 1e-7) try {

  auto norm = [&dot](arma::Col<coeff_t> const &v) {
    return std::sqrt(hydra::real(dot(v, v)));
  };

  auto tmatrix = Tmatrix();

  // Initialize Lanczos vectors and tmatrix
  arma::Col<coeff_t> v1, w;
  try {
    v1 = v0;
    w.resize(v0.size());
    w.zeros();
    v0.zeros();
  } catch (...) {
    throw(std::runtime_error("Cannot allocate Lanczos vectors"));
  }

  double alpha = 0.;
  double beta = 0.;

  // Normalize start vector or return if norm is zero
  coeff_t v1_norm = norm(v1);
  if (std::abs(v1_norm) > 1e-12) {
    v1 /= v1_norm;
  } else {
    return lanczos_result_t();
  }

  // Main Lanczos loop
  int64_t iteration = 0;
  std::string criterion;
  while (!converged(tmatrix)) {
    operation(v1);
    lanczos_step(v0, v1, w, alpha, beta, mult, dot);

    tmatrix.append(alpha, beta);
    tmatrix.print_log();

    ++iteration;

    // Finish if Lanczos sequence is exhausted
    if (std::abs(beta) > deflation_tol) {
      v1 /= beta;
    } else {
      criterion = "deflated";
      break;
    }

    if (iteration >= max_iterations) {
      criterion = "maxiterations";
      break;
    }
  }

  if (converged(tmatrix)) {
    criterion = "converged";
  }

  return lanczos_result_t{tmatrix.alphas(), tmatrix.betas(),
                          tmatrix.eigenvalues(), iteration, criterion};
} catch (...) {
  HydraRethrow("Error performing Lanczos algorithm");
  return lanczos_result_t();
}

} // namespace hydra::lanczos
