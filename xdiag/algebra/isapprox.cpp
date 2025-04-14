// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "isapprox.hpp"

namespace xdiag {
bool isapprox(double x, double y, double rtol, double atol) {
  return std::abs(x - y) <= (atol + rtol * std::abs(y));
}
bool isapprox(complex x, complex y, double rtol, double atol) {
  return std::abs(x - y) <= (atol + rtol * std::abs(y));
}

bool isapprox(arma::vec const &v, arma::vec const &w, double rtol,
              double atol) {
  return arma::approx_equal(v, w, "both", atol, rtol);
}
bool isapprox(arma::cx_vec const &v, arma::cx_vec const &w, double rtol,
              double atol) {
  return arma::approx_equal(v, w, "both", atol, rtol);
}

bool isapprox(arma::mat const &A, arma::mat const &B, double rtol,
              double atol) {
  return arma::approx_equal(A, B, "both", atol, rtol);
}
bool isapprox(arma::cx_mat const &A, arma::cx_mat const &B, double rtol,
              double atol) {
  return arma::approx_equal(A, B, "both", atol, rtol);
}

} // namespace xdiag
