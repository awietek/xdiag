#include "orthogonalize.hpp"

#include <xdiag/utils/logger.hpp>

namespace xdiag {

template <typename vec_t, typename mat_t>
void orthogonalize_inplace(vec_t &&v, mat_t const &Q, int64_t max_col,
                           int iterations) {
  if (max_col < 0) {
    max_col = Q.n_cols;
  } else if (max_col > (int64_t)Q.n_cols) {
    Log.err("Error in orthogonalize_inplace: max_col larger than number of "
            "matrix columns");
  }

  for (int iter = 0; iter < iterations; ++iter) {
    for (int64_t j = 0; j < max_col; ++j) {
      auto proj = cdot(Q.col(j), v);
      v -= proj * Q.col(j);
    }
  }
}

template void orthogonalize_inplace(arma::Col<double> &,
                                    arma::Mat<double> const &, int64_t, int);
template void orthogonalize_inplace(arma::Col<complex> &,
                                    arma::Mat<complex> const &, int64_t, int);
template void orthogonalize_inplace(arma::subview_col<double> &&,
                                    arma::Mat<double> const &, int64_t, int);
template void orthogonalize_inplace(arma::subview_col<complex> &&,
                                    arma::Mat<complex> const &, int64_t, int);

template <typename vec_t, typename mat_t>
vec_t orthogonalize(vec_t const &v, mat_t const &Q, int64_t max_col,
                    int iterations) {
  vec_t w = v;
  orthogonalize_inplace(w, Q, max_col, iterations);
  return w;
}

template arma::Col<double> orthogonalize(arma::Col<double> const &v,
                                         arma::Mat<double> const &Q,
                                         int64_t max_col, int iterations);
template arma::Col<complex> orthogonalize(arma::Col<complex> const &v,
                                          arma::Mat<complex> const &Q,
                                          int64_t max_col, int iterations);
template arma::subview_col<double>
orthogonalize(arma::subview_col<double> const &v, arma::Mat<double> const &Q,
              int64_t max_col, int iterations);
template arma::subview_col<complex>
orthogonalize(arma::subview_col<complex> const &v, arma::Mat<complex> const &Q,
              int64_t max_col, int iterations);

} // namespace xdiag
