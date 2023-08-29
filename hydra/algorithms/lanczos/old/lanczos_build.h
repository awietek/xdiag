#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/algorithms/lanczos/tmatrix.h>
#include <hydra/common.h>
#include <hydra/utils/logger.h>

namespace hydra {

// Lanczos implementation building a single vector
template <class coeff_t, class multiply_f>
std::pair<Tmatrix, arma::Mat<coeff_t>>
lanczos_build(multiply_f mult, arma::Col<coeff_t> &v0,
              arma::Col<coeff_t> coefficients, double deflation_tol = 1e-7) {

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
  int n_iterations = coefficients.size();

  arma::Col<coeff_t> vector(v0.size(), arma::fill::zeros);

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
    return {tmatrix, vector};
  }

  // Main Lanczos loop
  for (int iteration = 0; iteration < n_iterations; ++iteration) {

    // Build vector from linear combination
    vector += coefficients(iteration) * v1;

    lanczos_step(v0, v1, w, alpha, beta, mult, dot);
    tmatrix.append(alpha, beta);
    tmatrix.print_log();

    // Finish if Lanczos sequence is exhausted
    if (std::abs(beta) > deflation_tol) {
      v1 /= beta;
    } else {
      break;
    }
  }

  return {tmatrix, vector};
}

// Lanczos implementation building a single vector
template <class coeff_t, class multiply_f>
std::pair<Tmatrix, arma::Col<coeff_t>>
lanczos_build(multiply_f mult, arma::Col<coeff_t> &v0,
             arma::Col<real_t<coeff_t>> coefficients,
             double deflation_tol = 1e-7) {

  arma::Col<coeff_t> coefficients_cplx(coefficients.size());

  for (arma::uword i = 0; i < coefficients.size(); ++i) {
    coefficients_cplx(i) = (coeff_t)coefficients(i);
  }
  return lanczos_build(mult, v0, coefficients_cplx, deflation_tol);
}

// Lanczos implementation building multiple vectors
template <class coeff_t, class multiply_f>
std::pair<Tmatrix, arma::Mat<coeff_t>>
lanczos_build(multiply_f mult, arma::Col<coeff_t> &v0,
             arma::Mat<coeff_t> coefficients, double deflation_tol = 1e-7) {

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
  auto n_vectors = coefficients.n_cols;
  int n_iterations = coefficients.n_rows;

  arma::Mat<coeff_t> vectors(v0.size(), n_vectors, arma::fill::zeros);

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
    return {tmatrix, vectors};
  }

  // Main Lanczos loop
  for (int iteration = 0; iteration < n_iterations; ++iteration) {

    // Build vectors from linear combination
    for (arma::uword v_idx = 0; v_idx < n_vectors; ++v_idx) {
      vectors.col(v_idx) += coefficients(iteration, v_idx) * v1;
    }

    lanczos_step(v0, v1, w, alpha, beta, mult, dot);
    tmatrix.append(alpha, beta);
    tmatrix.print_log();

    // Finish if Lanczos sequence is exhausted
    if (std::abs(beta) > deflation_tol) {
      v1 /= beta;
    } else {
      break;
    }
  }

  return {tmatrix, vectors};
}

// Lanczos implementation building multiple vectors
template <class coeff_t, class multiply_f>
std::pair<Tmatrix, arma::Mat<coeff_t>>
lanczos_build(multiply_f mult, arma::Col<coeff_t> &v0,
              arma::Mat<real_t<coeff_t>> coefficients,
              double deflation_tol = 1e-7) {

  arma::Mat<coeff_t> coefficients_cplx(coefficients.n_rows,
                                       coefficients.n_cols);

  for (arma::uword j = 0; j < coefficients.n_cols; ++j) {
    for (arma::uword i = 0; i < coefficients.n_rows; ++i) {
      coefficients_cplx(i, j) = (coeff_t)coefficients(i, j);
    }
  }

  return lanczos_build(mult, v0, coefficients_cplx, deflation_tol);
}

} // namespace hydra
