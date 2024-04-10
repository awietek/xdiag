#pragma once

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/algorithms/lanczos/tmatrix.hpp>
#include <xdiag/algorithms/gram_schmidt/orthogonalize.hpp>
#include <xdiag/common.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/print_macro.hpp>

namespace xdiag {

template <typename coeff_t, class multiply_f, class dot_f>
inline void lanczos_step(arma::Col<coeff_t> &v0, arma::Col<coeff_t> &v1,
                         arma::Col<coeff_t> &w, double &alpha, double &beta,
                         multiply_f mult, dot_f dot) try {

  auto norm = [&dot](arma::Col<coeff_t> const &v) {
    return std::sqrt(xdiag::real(dot(v, v)));
  };
  
  mult(v1, w); // MVM
  alpha = xdiag::real(dot(v1, w)); 
  w -= alpha * v1;
  w -= beta * v0;
  v0 = v1;
  v1 = w;
  beta = norm(v1);
} catch (...) {
  XDiagRethrow("Error performing Lanczos step");
}

template <typename coeff_t, class multiply_f, class dot_f>
inline void lanczos_step_ortho(arma::Col<coeff_t> &v0, arma::Col<coeff_t> &v1,
			       arma::Col<coeff_t> &w, double &alpha, double &beta,
			       multiply_f mult, dot_f dot, arma::Mat<coeff_t> const& V,
			       int64_t iteration) try {

  auto norm = [&dot](arma::Col<coeff_t> const &v) {
    return std::sqrt(xdiag::real(dot(v, v)));
  };

  mult(v1, w); // MVM
  alpha = xdiag::real(dot(v1, w));
  orthogonalize_inplace(w, V, iteration);
  v0 = v1;
  v1 = w;
  beta = norm(v1);
} catch (...) {
  XDiagRethrow("Error performing Lanczos step with orthogonalization");
}

  

} // namespace xdiag
