#pragma once

#include <lila/all.h>

#include <hydra/linalg/lanczos/tmatrix.h>

namespace hydra {

template <class coeff_t, class multiply_f, class dot_f, class convergence_f>
std::pair<Tmatrix<coeff_t>, lila::Matrix<coeff_t>>
LanczosGeneric(multiply_f mult, lila::Vector<coeff_t> &v0, dot_f dot,
               convergence_f converged, int max_iterations = 1000,
               lila::Matrix<coeff_t> coefficients = lila::Matrix<coeff_t>(),
               double deflation_tol = 1e-7) {
  using namespace lila;
  using real = real_t<coeff_t>;

  auto norm = [&dot](lila::Vector<coeff_t> const & v) {
    return std::sqrt(lila::real(dot(v, v)));
  };

  auto tmatrix = Tmatrix<coeff_t>();
  auto n_vectors = coefficients.nrows();
  auto vectors = Zeros<coeff_t>(v0.size(), n_vectors);

  bool only_build_vectors = n_vectors > 0;
  if (only_build_vectors)
    max_iterations = coefficients.ncols();

  // Initialize Lanczos vectors and tmatrix
  auto v1 = v0;
  auto w = ZerosLike(v0);
  Zeros(v0);
  real alpha = 0.;
  real beta = 0.;

  // Normalize start vector or return if norm is zero
  coeff_t v1_norm = norm(v1);
  if (!lila::close(v1_norm, (coeff_t)0.)){
    Scale(1. / v1_norm, v1);
  }  else {
    return {tmatrix, vectors};
  }

  // Main Lanczos loop
  int iteration = 0;
  while (!converged(tmatrix) || only_build_vectors) {
    // Build vectors
    for (int v_idx = 0; v_idx < coefficients.nrows(); ++v_idx)
      Add(v1, vectors(ALL, v_idx), coefficients(v_idx, iteration));

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

} // namespace hydra
