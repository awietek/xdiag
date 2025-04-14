// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/algorithms/lanczos/lanczos.hpp>
#include <xdiag/algorithms/lanczos/lanczos_convergence.hpp>
#include <xdiag/algorithms/lanczos/tmatrix.hpp>

namespace xdiag {

struct exp_sym_v_result_t {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  int64_t niterations;
  std::string criterion;
};

  template <typename coeff_t, class multiply_f, class dot_f>
exp_sym_v_result_t
  exp_sym_v(multiply_f mult, dot_f dot, arma::Col<coeff_t> &X, coeff_t tau,
          double precision = 1e-12, double shift = 0, bool normalize = false,
          int64_t max_iterations = 1000, double deflation_tol = 1e-7) try {

  auto norm_f = [&dot](arma::Col<coeff_t> const &v) {
    return std::sqrt(xdiag::real(dot(v, v)));
  };
    
  double norm = norm_f(X);
  auto v0 = X;

  auto converged = [precision, tau, norm](Tmatrix const &tmat) {
    return lanczos::converged_time_evolution(tmat, tau, precision, norm);
  };

  auto operation_void = [](arma::Col<coeff_t> const &) {};
  
  auto r = lanczos::lanczos(mult, dot, converged, operation_void, v0,
                            max_iterations, deflation_tol);
  
  double e0 = r.eigenvalues(0);

  // Compute the tridiagonal matrix
  arma::mat tmat = arma::diagmat(r.alphas);
  if (r.alphas.n_rows > 1) {
    tmat += arma::diagmat(r.betas.head(r.betas.size() - 1), 1) +
            arma::diagmat(r.betas.head(r.betas.size() - 1), -1);
  }

  // Subtract shift from diagonal
  if (shift != 0.) {
    for (int64_t i = 0; i < (int64_t)tmat.n_cols; ++i) {
      tmat(i, i) -= shift;
    }
  }

  arma::Mat<coeff_t> texp = arma::expmat(tau * tmat);
  arma::Col<coeff_t> linear_combination = texp.col(0);

  v0 = X;
  X.zeros();
  int64_t iter = 0;
  auto mult2 = [&](arma::Col<coeff_t> const &v, arma::Col<coeff_t> &w) {
    mult(v, w);
    ++iter;
  };
  auto operation = [&linear_combination, &iter,
                    &X](arma::Col<coeff_t> const &v) {
    X += linear_combination(iter) * v;
  };

  lanczos::lanczos(mult2, dot, converged, operation, v0, max_iterations,
                   deflation_tol);

  if (!normalize) {
    X *= norm;
  } else {
    double nrm = norm_f(X);
    X /= nrm;
  }

  return {r.alphas, r.betas, r.eigenvalues, r.niterations, r.criterion};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
