#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/common.h>
#include <hydra/utils/precision.h>

namespace hydra {

inline bool close(double const &x, double const &y,
                  real_t<double> atol = atol<double>::val(),
                  real_t<double> rtol = rtol<double>::val()) {
  return (std::abs(x - y) <= (atol + rtol * std::abs(y)));
}

inline bool close(complex const &x, complex const &y,
                  real_t<complex> atol = atol<complex>::val(),
                  real_t<complex> rtol = rtol<complex>::val()) {
  return (std::abs(x - y) <= (atol + rtol * std::abs(y)));
}

template <class coeff_t>
inline bool close(arma::Col<coeff_t> const &vec1,
                  arma::Col<coeff_t> const &vec2,
                  real_t<coeff_t> atol = atol<coeff_t>::val(),
                  real_t<coeff_t> rtol = rtol<coeff_t>::val()) {
  assert(vec1.size() == vec2.size());
  return arma::approx_equal(vec1, vec2, "both", atol, rtol);
}

template <class coeff_t>
inline bool close(arma::Mat<coeff_t> const &mat1,
                  arma::Mat<coeff_t> const &mat2,
                  real_t<coeff_t> atol = atol<coeff_t>::val(),
                  real_t<coeff_t> rtol = rtol<coeff_t>::val()) {
  assert(mat1.n_rows == mat2.n_rows);
  assert(mat1.n_cols == mat2.n_cols);
  return arma::approx_equal(mat1, mat2, "both", atol, rtol);
}

} // namespace hydra
