//
// Created by Luke Staszewski on 07.02.23.
//
#include <algorithm>
#include <cmath>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/common.h>

namespace xdiag {

template <class coeff_t> static double log2abs(coeff_t x) {
  double value;
  if (x == 0.0)
    value = -1e30;
  else
    value = std::log(std::abs(x)) / std::log(2.0);

  return value;
}

template <typename coeff_t>
arma::Mat<coeff_t> expm(arma::Mat<coeff_t> const &A, coeff_t alpha = 1.) {
  using mat_t = arma::Mat<coeff_t>;
  arma::uword n = A.n_rows;
  assert(n == A.n_cols);

  const int q = 6;
  mat_t A2 = alpha * A;
  auto a_norm = arma::norm(A2, "inf");
  int ee = (int)log2abs(a_norm) + 1;
  int s = std::max(0, ee + 1);
  double t = 1.0 / pow(2.0, s);

  A2 = t * A2;
  mat_t X = A2;
  double c = 0.5;
  mat_t E = arma::eye<mat_t>(n, n);
  E += c * A2;
  mat_t D = arma::eye<mat_t>(n, n);
  D -= c * A2;

  for (int k = 2; k <= q; k++) {
    c *= (q - k + 1) / (double)(k * (2 * q - k + 1));
    X = X * A2;
    E += c * X;

    if (k % 2 == 0)
      D += c * X;
    else
      D -= c * X;
  }

  // inverse of D
  E = E * D.i();
  E = arma::powmat(E, (int)std::pow(2, s));

  return E;
}

template arma::Mat<double> expm(arma::Mat<double> const &A, double alpha = 1.);
template arma::Mat<complex> expm(arma::Mat<complex> const &A,
                                 complex alpha = 1.);
} // namespace xdiag
