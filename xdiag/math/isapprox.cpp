// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "isapprox.hpp"

#include <cmath>

namespace xdiag {

bool isapprox(double a, double b, double rtol, double atol) {
  return std::abs(a - b) <= atol + rtol * std::max(std::abs(a), std::abs(b));
}

bool isapprox(complex a, complex b, double rtol, double atol) {
  return std::abs(a - b) <= atol + rtol * std::max(std::abs(a), std::abs(b));
}

template <typename coeff_t>
bool isapprox(arma::Col<coeff_t> const &v, arma::Col<coeff_t> const &w,
              double rtol, double atol) {
  return arma::approx_equal(v, w, "both", atol, rtol);
}

template <typename coeff_t>
bool isapprox(arma::Mat<coeff_t> const &v, arma::Mat<coeff_t> const &w,
              double rtol, double atol) {
  return arma::approx_equal(v, w, "both", atol, rtol);
}

template bool isapprox(arma::vec const &, arma::vec const &, double, double);
template bool isapprox(arma::cx_vec const &, arma::cx_vec const &, double,
                       double);
template bool isapprox(arma::mat const &, arma::mat const &, double, double);
template bool isapprox(arma::cx_mat const &, arma::cx_mat const &, double,
                       double);
} // namespace xdiag
