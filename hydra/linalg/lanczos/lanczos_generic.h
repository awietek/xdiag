#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/common.h>
#include <hydra/linalg/lanczos/tmatrix.h>
#include <hydra/utils/logger.h>

namespace hydra {

// Generic Lanczos implementation building multiple vectors
template <class coeff_t, class multiply_f, class dot_f, class convergence_f>
std::pair<Tmatrix, arma::Mat<coeff_t>>
LanczosGeneric(multiply_f mult, arma::Col<coeff_t> &v0, dot_f dot,
               convergence_f converged,
               arma::Mat<coeff_t> coefficients = arma::Mat<coeff_t>(),
               int max_iterations = 1000, double deflation_tol = 1e-7) {

  using real = real_t<coeff_t>;

  auto norm = [&dot](arma::Col<coeff_t> const &v) {
    return std::sqrt(hydra::real(dot(v, v)));
  };

  auto tmatrix = Tmatrix();
  auto n_vectors = coefficients.n_cols;

  arma::Mat<coeff_t> vectors;
  if (n_vectors > 0)
    vectors = arma::zeros<arma::Mat<coeff_t>>(v0.size(), n_vectors);

  bool only_build_vectors = n_vectors > 0;
  if (only_build_vectors)
    max_iterations = coefficients.n_rows;

  // Initialize Lanczos vectors and tmatrix
  auto v1 = v0;
  arma::Col<coeff_t> w(v0.size(), arma::fill::zeros);
  v0.zeros();
  real alpha = 0.;
  real beta = 0.;

  // Normalize start vector or return if norm is zero
  coeff_t v1_norm = norm(v1);
  if (std::abs(v1_norm) > 1e-12) {
    v1 /= v1_norm;
  } else {
    return {tmatrix, vectors};
  }

  // Main Lanczos loop
  int iteration = 0;
  while (!converged(tmatrix) || only_build_vectors) {
    // Build vectors from linear combination
    for (arma::uword v_idx = 0; v_idx < coefficients.n_cols; ++v_idx) {
      vectors.col(v_idx) += coefficients(iteration, v_idx) * v1;
    }

    // Lanczos recursion
    mult(v1, w); // MVM
    alpha = hydra::real(dot(v1, w));
    w -= alpha * v1;
    w -= beta * v0;
    v0 = v1;
    v1 = w;
    beta = norm(v1);
    tmatrix.append(alpha, beta);

    auto eigs = tmatrix.eigenvalues();
    Log(2, "alpha: {:.16f}", alpha);
    Log(2, "beta: {:.16f}", beta);
    if (eigs.size() == 1) {
      Log(2, "eigs: {:.16f}", eigs(0));
    } else if (eigs.size() == 2) {
      Log(2, "eigs: {:.16f} {:.16f}", eigs(0), eigs(1));
    } else {
      Log(2, "eigs: {:.16f} {:.16f} {:.16f}", eigs(0), eigs(1), eigs(2));
    }

    // Finish if Lanczos sequence is exhausted
    if (std::abs(beta) > deflation_tol) {
      v1 /= beta;
    } else
      break;

    ++iteration;
    if (iteration >= max_iterations)
      break;
  }

  return {tmatrix, vectors};
}

// Generic Lanczos implementation building multiple vectors (real coefficients)
template <class multiply_f, class dot_f, class convergence_f>
std::pair<Tmatrix, arma::Mat<complex>>
LanczosGeneric(multiply_f mult, arma::Col<complex> &v0, dot_f dot,
               convergence_f converged,
               arma::Mat<double> coefficients = arma::Mat<double>(),
               int max_iterations = 1000, double deflation_tol = 1e-7) {
  return LanczosGeneric(mult, v0, dot, converged, coefficients, max_iterations,
                        deflation_tol);
}

// Generic Lanczos implementation building one vector (real coefficients)
template <class coeff_t, class multiply_f, class dot_f, class convergence_f>
std::pair<Tmatrix, arma::Col<coeff_t>>
LanczosGeneric(multiply_f mult, arma::Col<coeff_t> &v0, dot_f dot,
               convergence_f converged,
               arma::Col<real_t<coeff_t>> const &coefficients,
               int max_iterations = 1000, double deflation_tol = 1e-7) {
  if (coefficients.size() > 0) {
    auto coefficients_mat = arma::Mat<coeff_t>(Complex(coefficients));
    auto [tmat, vecs] =
        LanczosGeneric(mult, v0, dot, converged, coefficients_mat,
                       max_iterations, deflation_tol);
    return {tmat, vecs.col(0)};
  } else {
    auto [tmat, vecs] =
        LanczosGeneric(mult, v0, dot, converged, arma::Mat<coeff_t>(),
                       max_iterations, deflation_tol);
    return {tmat, vecs.col(0)};
  }
}

// Generic Lanczos implementation building one vector (complex coefficients)
template <class coeff_t, class multiply_f, class dot_f, class convergence_f>
std::pair<Tmatrix, arma::Col<coeff_t>>
LanczosGeneric(multiply_f mult, arma::Col<coeff_t> &v0, dot_f dot,
               convergence_f converged, arma::Col<coeff_t> const &coefficients,
               int max_iterations = 1000, double deflation_tol = 1e-7) {
  if (coefficients.size() > 0) {
    auto coefficients_mat = arma::Mat<coeff_t>(coefficients);
    auto [tmat, vecs] =
        LanczosGeneric(mult, v0, dot, converged, coefficients_mat,
                       max_iterations, deflation_tol);
    return {tmat, vecs.col(0)};
  } else {
    auto [tmat, vecs] =
        LanczosGeneric(mult, v0, dot, converged, arma::Mat<coeff_t>(),
                       max_iterations, deflation_tol);
    return {tmat, vecs.col(0)};
  }
}

} // namespace hydra
