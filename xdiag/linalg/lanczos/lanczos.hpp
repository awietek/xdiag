// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/linalg/lanczos/lanczos_step.hpp>
#include <xdiag/linalg/lanczos/tmatrix.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::lanczos {

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
  // Below this norm the start vector is treated as zero (unusable). This is a
  // separate notion from deflation_tol, which detects an exhausted sequence.
  constexpr double start_vector_tol = 1e-12;
  auto tmatrix = Tmatrix();

  // Initialize Lanczos vectors and tmatrix
  arma::Col<coeff_t> v1, w;
  try {
    v1 = v0;
    w.resize(v0.size());
    w.zeros();
    v0.zeros();
  } catch (...) {
    XDIAG_THROW("Cannot allocate Lanczos vectors");
  }

  double alpha = 0.;
  double beta = 0.;

  // Normalize start vector or return if norm is zero
  double v1_norm = lanczos_norm(v1, dot);
  if (v1_norm > start_vector_tol) {
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
}
XDIAG_CATCH

} // namespace xdiag::lanczos
