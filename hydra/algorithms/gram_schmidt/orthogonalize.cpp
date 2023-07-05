#include "orthogonalize.h"

namespace hydra {
template <typename coeff_t>
void orthogonalize_inplace(arma::Col<coeff_t> &v, arma::Mat<coeff_t> const &Q,
                           int64_t max_col, int iterations) {
  if (max_col < 0) {
    max_col = Q.n_cols;
  }
  for (int iter = 0; iter < iterations; ++iter) {
    for (int j = 0; j < max_col; ++j) {
      coeff_t proj = cdot(Q.col(j), v);
      v -= proj * Q.col(j);
    }
  }
}

template void orthogonalize_inplace(arma::Col<double> &,
                                    arma::Mat<double> const &, int64_t, int);
template void orthogonalize_inplace(arma::Col<complex> &,
                                    arma::Mat<complex> const &, int64_t, int);

template <typename coeff_t>
arma::Col<coeff_t> orthogonalize(arma::Col<coeff_t> const &v,
                                 arma::Mat<coeff_t> const &Q, int64_t max_col,
                                 int iterations) {
  arma::Col<coeff_t> w = v;
  orthogonalize_inplace(w, Q, max_col, iterations);
  return w;
}

template arma::Col<double> orthogonalize(arma::Col<double> const &v,
                                         arma::Mat<double> const &Q,
                                         int64_t max_col, int iterations);
template arma::Col<complex> orthogonalize(arma::Col<complex> const &v,
                                          arma::Mat<complex> const &Q,
                                          int64_t max_col, int iterations);
} // namespace hydra
