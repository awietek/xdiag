// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/linalg/gram_schmidt/orthogonalize.hpp>
#include <xdiag/linalg/lanczos/tmatrix.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

// Norm induced by the (possibly distributed) dot product.
template <typename coeff_t, class dot_f>
inline double lanczos_norm(arma::Col<coeff_t> const &v, dot_f dot) {
  return std::sqrt(xdiag::real(dot(v, v)));
}

template <typename coeff_t, class multiply_f, class dot_f>
inline void lanczos_step(arma::Col<coeff_t> &v0, arma::Col<coeff_t> &v1,
                         arma::Col<coeff_t> &w, double &alpha, double &beta,
                         multiply_f mult, dot_f dot) try {
  mult(v1, w); // MVM
  alpha = xdiag::real(dot(v1, w));
  w -= alpha * v1;
  w -= beta * v0;
  // Rotate the three-term recurrence. v0 may alias external (State) memory and
  // must be assigned element-wise; v1 and w are locally owned so we swap them
  // to avoid a second full-length copy.
  v0 = v1;
  v1.swap(w);
  beta = lanczos_norm(v1, dot);
}
XDIAG_CATCH

template <typename coeff_t, class multiply_f, class dot_f>
inline void lanczos_step_ortho(arma::Col<coeff_t> &v0, arma::Col<coeff_t> &v1,
                               arma::Col<coeff_t> &w, double &alpha,
                               double &beta, multiply_f mult, dot_f dot,
                               arma::Mat<coeff_t> const &V,
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
}
XDIAG_CATCH

} // namespace xdiag
