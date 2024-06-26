#pragma once

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/common.hpp>
#include <xdiag/utils/precision.hpp>

namespace xdiag {

inline bool close(double const &x, double const &y,
                  real_t<double> atol = numeric::atol<double>::val(),
                  real_t<double> rtol = numeric::rtol<double>::val()) {
  return (std::abs(x - y) <= (atol + rtol * std::abs(y)));
}

inline bool close(complex const &x, complex const &y,
                  real_t<complex> atol = numeric::atol<complex>::val(),
                  real_t<complex> rtol = numeric::rtol<complex>::val()) {
  return (std::abs(x - y) <= (atol + rtol * std::abs(y)));
}

template <class coeff_t>
inline bool close(arma::Col<coeff_t> const &vec1,
                  arma::Col<coeff_t> const &vec2,
                  real_t<coeff_t> atol = numeric::atol<coeff_t>::val(),
                  real_t<coeff_t> rtol = numeric::rtol<coeff_t>::val()) {
  assert(vec1.size() == vec2.size());
  return arma::approx_equal(vec1, vec2, "both", atol, rtol);
}

template <class coeff_t>
inline bool close(arma::Mat<coeff_t> const &mat1,
                  arma::Mat<coeff_t> const &mat2,
                  real_t<coeff_t> atol = numeric::atol<coeff_t>::val(),
                  real_t<coeff_t> rtol = numeric::rtol<coeff_t>::val()) {
  assert(mat1.n_rows == mat2.n_rows);
  assert(mat1.n_cols == mat2.n_cols);
  return arma::approx_equal(mat1, mat2, "both", atol, rtol);
}

} // namespace xdiag
