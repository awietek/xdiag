#pragma once

#include <lila/all.h>
#include <hydra/linalg/lanczos/tmatrix.h>
#include <hydra/common.h>

namespace hydra {

// Generic Lanczos implementation building multiple vectors
template <class coeff_t, class multiply_f, class dot_f, class convergence_f>
std::pair<Tmatrix, lila::Matrix<coeff_t>>
LanczosGeneric(multiply_f mult, lila::Vector<coeff_t> &v0, dot_f dot,
               convergence_f converged,
               lila::Matrix<coeff_t> coefficients = lila::Matrix<coeff_t>(),
               int max_iterations = 1000, double deflation_tol = 1e-7) {

  using namespace lila;
  using real = real_t<coeff_t>;

  auto norm = [&dot](lila::Vector<coeff_t> const &v) {
    return std::sqrt(lila::real(dot(v, v)));
  };

  auto tmatrix = Tmatrix();
  auto n_vectors = coefficients.n();

  lila::Matrix<coeff_t> vectors;
  if (n_vectors > 0)
    vectors = Zeros<coeff_t>(v0.size(), n_vectors);

  bool only_build_vectors = n_vectors > 0;
  if (only_build_vectors)
    max_iterations = coefficients.m();

  // Initialize Lanczos vectors and tmatrix
  auto v1 = v0;
  auto w = ZerosLike(v0);
  Zeros(v0);
  real alpha = 0.;
  real beta = 0.;

  // Normalize start vector or return if norm is zero
  coeff_t v1_norm = norm(v1);
  if (!lila::close(v1_norm, (coeff_t)0.)) {
    Scale(1. / v1_norm, v1);
  } else {
    return {tmatrix, vectors};
  }

  // Main Lanczos loop
  int iteration = 0;
  while (!converged(tmatrix) || only_build_vectors) {
    // Build vectors from linear combination
    for (int v_idx = 0; v_idx < coefficients.n(); ++v_idx) {
      Add(v1, vectors({0, vectors.m()}, v_idx), coefficients(iteration, v_idx));
    }

    // Lanczos recursion
    mult(v1, w); // MVM
    alpha = lila::real(dot(v1, w));
    Add(v1, w, -(coeff_t)alpha); // w -= alpha*v1;
    Add(v0, w, -(coeff_t)beta);  // w -= beta*v0;
    v0 = v1;
    v1 = w;
    beta = norm(v1);
    tmatrix.append(alpha, beta);

    // Finish if Lanczos sequence is exhausted
    if (std::abs(beta) > deflation_tol) {
      Scale(1. / (coeff_t)beta, v1);
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
std::pair<Tmatrix, lila::Matrix<complex>>
LanczosGeneric(multiply_f mult, lila::Vector<complex> &v0, dot_f dot,
               convergence_f converged,
               lila::Matrix<double> coefficients = lila::Matrix<double>(),
               int max_iterations = 1000, double deflation_tol = 1e-7) {
  return LanczosGeneric(mult, v0, dot, converged, Complex(coefficients),
                        max_iterations, deflation_tol);
}

// Generic Lanczos implementation building one vector (real coefficients)
template <class coeff_t, class multiply_f, class dot_f, class convergence_f>
std::pair<Tmatrix, lila::Vector<coeff_t>>
LanczosGeneric(multiply_f mult, lila::Vector<coeff_t> &v0, dot_f dot,
               convergence_f converged,
               lila::Vector<lila::real_t<coeff_t>> const &coefficients,
               int max_iterations = 1000, double deflation_tol = 1e-7) {
  if (coefficients.size() > 0) {
    auto coefficients_mat = lila::Matrix<coeff_t>(Complex(coefficients));
    auto [tmat, vecs] =
        LanczosGeneric(mult, v0, dot, converged, coefficients_mat,
                       max_iterations, deflation_tol);
    return {tmat, vecs.col(0)};
  } else {
    auto [tmat, vecs] =
        LanczosGeneric(mult, v0, dot, converged, lila::Matrix<coeff_t>(),
                       max_iterations, deflation_tol);
    return {tmat, vecs.col(0)};
  }
}

// Generic Lanczos implementation building one vector (complex coefficients)
template <class coeff_t, class multiply_f, class dot_f, class convergence_f>
std::pair<Tmatrix, lila::Vector<coeff_t>>
LanczosGeneric(multiply_f mult, lila::Vector<coeff_t> &v0, dot_f dot,
               convergence_f converged,
               lila::Vector<coeff_t> const &coefficients,
               int max_iterations = 1000, double deflation_tol = 1e-7) {
  if (coefficients.size() > 0) {
    auto coefficients_mat = lila::Matrix<coeff_t>(coefficients);
    auto [tmat, vecs] =
        LanczosGeneric(mult, v0, dot, converged, coefficients_mat,
                       max_iterations, deflation_tol);
    return {tmat, vecs.col(0)};
  } else {
    auto [tmat, vecs] =
        LanczosGeneric(mult, v0, dot, converged, lila::Matrix<coeff_t>(),
                       max_iterations, deflation_tol);
    return {tmat, vecs.col(0)};
  }
}

} // namespace hydra
