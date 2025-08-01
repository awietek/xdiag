// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/utils/logger.hpp>

namespace xdiag {

template <typename idx_t, typename coeff_t>
arma::Col<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                         arma::Col<coeff_t> const &vec_in) try {
  int64_t m = spmat.nrows;
  auto vec_out = arma::Col<coeff_t>(m);
  apply(spmat, vec_in, vec_out);
  return vec_out;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::vec apply(CSRMatrix<int32_t, double> const &spmat,
                         arma::vec const &vec_in);
template arma::vec apply(CSRMatrix<int64_t, double> const &spmat,
                         arma::vec const &vec_in);
template arma::cx_vec apply(CSRMatrix<int32_t, complex> const &spmat,
                            arma::cx_vec const &vec_in);
template arma::cx_vec apply(CSRMatrix<int64_t, complex> const &spmat,
                            arma::cx_vec const &vec_in);

template <typename idx_t, typename coeff_t>
arma::Mat<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                         arma::Mat<coeff_t> const &mat_in) try {
  int64_t m = spmat.nrows;
  int64_t k = mat_in.n_cols;
  auto mat_out = arma::Mat<coeff_t>(m, k);
  apply(spmat, mat_in, mat_out);
  return mat_out;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::mat apply(CSRMatrix<int32_t, double> const &spmat,
                         arma::mat const &mat_in);
template arma::mat apply(CSRMatrix<int64_t, double> const &spmat,
                         arma::mat const &mat_in);
template arma::cx_mat apply(CSRMatrix<int32_t, complex> const &spmat,
                            arma::cx_mat const &mat_in);
template arma::cx_mat apply(CSRMatrix<int64_t, complex> const &spmat,
                            arma::cx_mat const &mat_in);

template <typename idx_t, typename coeff_t>
XDIAG_API void apply(CSRMatrix<idx_t, coeff_t> const &A,
                     arma::Col<coeff_t> const &vec_in,
                     arma::Col<coeff_t> &vec_out) try {
  idx_t i0 = A.i0;
  if (!((i0 == 0) || (i0 == 1))) {
    XDIAG_THROW(fmt::format(
        "Invalid zero index i0. Must be either 0 or 1, but got i0={}", i0));
  }

  idx_t m = A.nrows;
  idx_t n = A.ncols;

  if ((m != vec_out.n_rows) || (n != vec_in.n_rows)) {
    XDIAG_THROW(fmt::format(
        "Incompatible sparse matrix and matrix dimensions A*x=y, sparse matrix "
        "A: ({} x {}), x: {}, y: {}",
        m, n, vec_in.n_rows, vec_out.n_rows))
  }

  if (i0 == 0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (idx_t i = 0; i < m; i++) {
      coeff_t zi = 0.0;
      idx_t j;
      for (idx_t ckey = A.rowptr[i]; ckey < A.rowptr[i + 1]; ckey++) {
        zi += A.data[ckey] * vec_in[A.col[ckey]];
      }
      vec_out[i] = zi;
    }
  } else { //  (i0 == 1)
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (idx_t i = 0; i < m; i++) {
      coeff_t zi = 0.0;
      idx_t j;
      for (idx_t idx = A.rowptr[i]; idx < A.rowptr[i + 1]; idx++) {
        zi += A.data[idx - 1] * vec_in[A.col[idx - 1] - 1];
      }
      vec_out[i] = zi;
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply(CSRMatrix<int32_t, double> const &, arma::vec const &,
                    arma::vec &);
template void apply(CSRMatrix<int64_t, double> const &, arma::vec const &,
                    arma::vec &);
template void apply(CSRMatrix<int32_t, complex> const &, arma::cx_vec const &,
                    arma::cx_vec &);
template void apply(CSRMatrix<int64_t, complex> const &, arma::cx_vec const &,
                    arma::cx_vec &);

template <typename idx_t, typename coeff_t>
XDIAG_API void apply(CSRMatrix<idx_t, coeff_t> const &A,
                     arma::Mat<coeff_t> const &mat_in,
                     arma::Mat<coeff_t> &mat_out) try {
  idx_t i0 = A.i0;
  if (!((i0 == 0) || (i0 == 1))) {
    XDIAG_THROW(fmt::format(
        "Invalid zero index i0. Must be either 0 or 1, but got i0={}", i0));
  }

  idx_t m = A.nrows;
  idx_t n = A.ncols;

  if ((m != mat_out.n_rows) || (n != mat_in.n_rows) ||
      (mat_in.n_cols != mat_out.n_cols)) {
    XDIAG_THROW(fmt::format(
        "Incompatible sparse matrix and matrix dimensions A*X=Y, sparse matrix "
        "A: ({} x {}), X: ({} x {}), Y: ({} x {})",
        m, n, mat_in.n_rows, mat_in.n_cols, mat_out.n_rows, mat_out.n_cols));
  }

  // implementation adapted from
  // https://stackoverflow.com/questions/76905042/what-is-wrong-with-my-sparse-matrix-multiple-vectors-spmm-product-function-for
  idx_t nvec = mat_out.n_cols;
  coeff_t *X = const_cast<coeff_t *>(mat_in.memptr());
  coeff_t *Y = mat_out.memptr();

  std::vector<coeff_t *> x(nvec);
  for (idx_t vec = 0; vec < nvec; vec++) {
    x[vec] = X + vec * n;
  }

  if (i0 == 0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (idx_t i = 0; i < m; i++) { // loop over rows
      coeff_t *y = Y + i;
      for (idx_t vec = 0; vec < nvec; vec++) {
        y[vec * m] = 0.;
      }
      for (idx_t jj = A.rowptr[i]; jj < A.rowptr[i + 1]; jj++) {
        coeff_t val = A.data[jj];
        idx_t j = A.col[jj];
        for (idx_t vec = 0; vec < nvec; vec++) {
          y[vec * m] += val * x[vec][j];
        }
      }
    }
  } else { // i0 == 1
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (idx_t i = 0; i < m; i++) { // loop over rows
      coeff_t *y = Y + i;
      for (idx_t vec = 0; vec < nvec; vec++) {
        y[vec * m] = 0.;
      }
      for (idx_t jj = A.rowptr[i]; jj < A.rowptr[i + 1]; jj++) {
        coeff_t val = A.data[jj - 1];
        idx_t j = A.col[jj - 1] - 1;
        for (idx_t vec = 0; vec < nvec; vec++) {
          y[vec * m] += val * x[vec][j];
        }
      }
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply(CSRMatrix<int32_t, double> const &, arma::mat const &,
                    arma::mat &);
template void apply(CSRMatrix<int64_t, double> const &, arma::mat const &,
                    arma::mat &);
template void apply(CSRMatrix<int32_t, complex> const &, arma::cx_mat const &,
                    arma::cx_mat &);
template void apply(CSRMatrix<int64_t, complex> const &, arma::cx_mat const &,
                    arma::cx_mat &);

} // namespace xdiag
