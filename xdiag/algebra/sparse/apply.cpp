// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef XDIAG_USE_SPARSE_MKL
#include <type_traits>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include <mkl_spblas.h>
#endif

namespace xdiag {

template <typename idx_t, typename coeff_csr_t, typename coeff_vec_t>
static arma::Col<coeff_vec_t> apply(CSRMatrix<idx_t, coeff_csr_t> const &spmat,
                                    arma::Col<coeff_vec_t> const &vec_in) try {
  static_assert(!(isreal<coeff_vec_t>() && !isreal<coeff_csr_t>()));
  int64_t m = spmat.nrows;
  auto vec_out = arma::Col<coeff_vec_t>(m);
  apply(spmat, vec_in, vec_out);
  return vec_out;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
arma::Col<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                         arma::Col<coeff_t> const &vec_in) try {
  return apply<idx_t, coeff_t, coeff_t>(spmat, vec_in);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template arma::vec apply(CSRMatrix<int32_t, double> const &, arma::vec const &);
template arma::vec apply(CSRMatrix<int64_t, double> const &, arma::vec const &);
template arma::cx_vec apply(CSRMatrix<int32_t, complex> const &,
                            arma::cx_vec const &);
template arma::cx_vec apply(CSRMatrix<int64_t, complex> const &,
                            arma::cx_vec const &);

template <typename idx_t>
arma::Col<complex> apply(CSRMatrix<idx_t, double> const &spmat,
                         arma::Col<complex> const &vec_in) try {
  return apply<idx_t, double, complex>(spmat, vec_in);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::cx_vec apply(CSRMatrix<int32_t, double> const &spmat,
                            arma::cx_vec const &vec_in);
template arma::cx_vec apply(CSRMatrix<int64_t, double> const &spmat,
                            arma::cx_vec const &vec_in);

template <typename idx_t, typename coeff_csr_t, typename coeff_mat_t>
static arma::Mat<coeff_mat_t> apply(CSRMatrix<idx_t, coeff_csr_t> const &spmat,
                                    arma::Mat<coeff_mat_t> const &mat_in) try {
  static_assert(!(isreal<coeff_mat_t>() && !isreal<coeff_csr_t>()));
  int64_t m = spmat.nrows;
  int64_t k = mat_in.n_cols;
  auto mat_out = arma::Mat<coeff_mat_t>(m, k);
  apply(spmat, mat_in, mat_out);
  return mat_out;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <typename idx_t, typename coeff_t>
arma::Mat<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                         arma::Mat<coeff_t> const &mat_in) try {
  return apply<idx_t, coeff_t, coeff_t>(spmat, mat_in);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template arma::mat apply(CSRMatrix<int32_t, double> const &, arma::mat const &);
template arma::mat apply(CSRMatrix<int64_t, double> const &, arma::mat const &);
template arma::cx_mat apply(CSRMatrix<int32_t, complex> const &,
                            arma::cx_mat const &);
template arma::cx_mat apply(CSRMatrix<int64_t, complex> const &,
                            arma::cx_mat const &);

template <typename idx_t>
arma::Mat<complex> apply(CSRMatrix<idx_t, double> const &spmat,
                         arma::Mat<complex> const &mat_in) try {
  return apply<idx_t, double, complex>(spmat, mat_in);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template arma::cx_mat apply(CSRMatrix<int32_t, double> const &spmat,
                            arma::cx_mat const &mat_in);
template arma::cx_mat apply(CSRMatrix<int64_t, double> const &spmat,
                            arma::cx_mat const &mat_in);

#ifdef XDIAG_USE_SPARSE_MKL
template <typename coeff_t>
static void apply_mkl(CSRMatrix<int64_t, coeff_t> const &A,
                      arma::Col<coeff_t> const &vec_in,
                      arma::Col<coeff_t> &vec_out) try {

  if ((A.nrows != vec_out.n_rows) || (A.ncols != vec_in.n_rows)) {
    XDIAG_THROW(fmt::format(
        "Incompatible sparse matrix and matrix dimensions A*x=y, sparse matrix "
        "A: ({} x {}), x: {}, y: {}",
        A.nrows, A.ncols, vec_in.n_rows, vec_out.n_rows))
  }

  sparse_matrix_t A_handle;
  sparse_status_t status;
  sparse_index_base_t i0 =
      (A.i0 == 0) ? SPARSE_INDEX_BASE_ZERO : SPARSE_INDEX_BASE_ONE;

  matrix_descr descr;
  descr.mode = SPARSE_FILL_MODE_LOWER;
  descr.diag = SPARSE_DIAG_NON_UNIT;

  MKL_INT nrows = A.nrows;
  MKL_INT ncols = A.ncols;
  MKL_INT *rowptr =
      reinterpret_cast<MKL_INT *>(const_cast<int64_t *>(A.rowptr.memptr()));
  MKL_INT *col =
      reinterpret_cast<MKL_INT *>(const_cast<int64_t *>(A.col.memptr()));

  // Real routines (double)
  if constexpr (std::is_same<coeff_t, double>::value) {
    double *data = const_cast<double *>(A.data.memptr());

    descr.type = A.ishermitian ? SPARSE_MATRIX_TYPE_SYMMETRIC
                               : SPARSE_MATRIX_TYPE_GENERAL;
    status = mkl_sparse_d_create_csr(&A_handle, i0, nrows, ncols, rowptr,
                                     rowptr + 1, col, data);
    if (status == SPARSE_STATUS_NOT_INITIALIZED) {
      // empty CSR matrix give -> return 0 vector
      vec_out.zeros();
      return;
    } else if (status != SPARSE_STATUS_SUCCESS) {
      XDIAG_THROW(
          fmt::format("Error creating MKL sparse matrix interface using "
                      "mkl_sparse_d_create_csr. MKL error code: {}",
                      (int)status));
    }
    mkl_sparse_set_mv_hint(A_handle, SPARSE_OPERATION_NON_TRANSPOSE, descr,
                           1000);
    status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, A_handle,
                             descr, vec_in.memptr(), 0.0, vec_out.memptr());
    if (status != SPARSE_STATUS_SUCCESS) {
      XDIAG_THROW(fmt::format("Error calling sparse matrix-vector multiply "
                              "mkl_sparse_d_mv (double). MKL error code: {}",
                              (int)status));
    }
  } else { // Complex routines
    complex *data = const_cast<complex *>(A.data.memptr());

    descr.type = A.ishermitian ? SPARSE_MATRIX_TYPE_HERMITIAN
                               : SPARSE_MATRIX_TYPE_GENERAL;

    status = mkl_sparse_z_create_csr(&A_handle, i0, nrows, ncols, rowptr,
                                     rowptr + 1, col, data);

    if (status == SPARSE_STATUS_NOT_INITIALIZED) {
      // empty CSR matrix give -> return 0 vector
      vec_out.zeros();
      return;
    } else if (status != SPARSE_STATUS_SUCCESS) {
      XDIAG_THROW(
          fmt::format("Error creating MKL sparse matrix interface using "
                      "mkl_sparse_z_create_csr. MKL error code: {}",
                      (int)status));
    }
    mkl_sparse_set_mv_hint(A_handle, SPARSE_OPERATION_NON_TRANSPOSE, descr,
                           1000);
    status = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, A_handle,
                             descr, vec_in.memptr(), 0.0, vec_out.memptr());

    if (status != SPARSE_STATUS_SUCCESS) {
      XDIAG_THROW(fmt::format("Error calling sparse matrix-vector multiply "
                              "mkl_sparse_z_mv (complex). MKL error code: {}",
                              (int)status));
    }
  }

  mkl_sparse_destroy(A_handle);
}
XDIAG_CATCH

template void apply_mkl(CSRMatrix<int64_t, double> const &, arma::vec const &,
                        arma::vec &);
template void apply_mkl(CSRMatrix<int64_t, complex> const &,
                        arma::cx_vec const &, arma::cx_vec &);
#endif

template <typename idx_t, typename coeff_csr_t, typename coeff_vec_t>
static void apply(CSRMatrix<idx_t, coeff_csr_t> const &A,
                  arma::Col<coeff_vec_t> const &vec_in,
                  arma::Col<coeff_vec_t> &vec_out) try {

  static_assert(!(isreal<coeff_vec_t>() && !isreal<coeff_csr_t>()));

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
      coeff_vec_t zi = 0.0;
      idx_t j;
      for (idx_t ckey = A.rowptr[i]; ckey < A.rowptr[i + 1]; ckey++) {
        zi += (coeff_vec_t)A.data[ckey] * vec_in[A.col[ckey]];
      }
      vec_out[i] = zi;
    }
  } else { //  (i0 == 1)
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (idx_t i = 0; i < m; i++) {
      coeff_vec_t zi = 0.0;
      idx_t j;
      for (idx_t idx = A.rowptr[i]; idx < A.rowptr[i + 1]; idx++) {
        zi += (coeff_vec_t)A.data[idx - 1] * vec_in[A.col[idx - 1] - 1];
      }
      vec_out[i] = zi;
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename coeff_t>
void apply(CSRMatrix<idx_t, coeff_t> const &spmat,
           arma::Col<coeff_t> const &vec_in, arma::Col<coeff_t> &vec_out) try {
#ifdef XDIAG_USE_SPARSE_MKL
  if constexpr (std::is_same<idx_t, int64_t>::value) {
    apply_mkl<coeff_t>(spmat, vec_in, vec_out);
  } else {
    apply<idx_t, coeff_t, coeff_t>(spmat, vec_in, vec_out);
  }
#else
  apply<idx_t, coeff_t, coeff_t>(spmat, vec_in, vec_out);
#endif
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

template <typename idx_t>
XDIAG_API void apply(CSRMatrix<idx_t, double> const &spmat,
                     arma::Col<complex> const &vec_in,
                     arma::Col<complex> &vec_out) try {
  apply<idx_t, double, complex>(spmat, vec_in, vec_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply(CSRMatrix<int32_t, double> const &, arma::cx_vec const &,
                    arma::cx_vec &);
template void apply(CSRMatrix<int64_t, double> const &, arma::cx_vec const &,
                    arma::cx_vec &);

template <typename idx_t, typename coeff_csr_t, typename coeff_mat_t>
static void apply(CSRMatrix<idx_t, coeff_csr_t> const &A,
                  arma::Mat<coeff_mat_t> const &mat_in,
                  arma::Mat<coeff_mat_t> &mat_out) try {
  static_assert(!(isreal<coeff_mat_t>() && !isreal<coeff_csr_t>()));

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
  coeff_mat_t *X = const_cast<coeff_mat_t *>(mat_in.memptr());
  coeff_mat_t *Y = mat_out.memptr();

  std::vector<coeff_mat_t *> x(nvec);
  for (idx_t vec = 0; vec < nvec; vec++) {
    x[vec] = X + vec * n;
  }

  if (i0 == 0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (idx_t i = 0; i < m; i++) { // loop over rows
      coeff_mat_t *y = Y + i;
      for (idx_t vec = 0; vec < nvec; vec++) {
        y[vec * m] = 0.;
      }
      for (idx_t jj = A.rowptr[i]; jj < A.rowptr[i + 1]; jj++) {
        coeff_mat_t val = (coeff_mat_t)A.data[jj];
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
      coeff_mat_t *y = Y + i;
      for (idx_t vec = 0; vec < nvec; vec++) {
        y[vec * m] = 0.;
      }
      for (idx_t jj = A.rowptr[i]; jj < A.rowptr[i + 1]; jj++) {
        coeff_mat_t val = (coeff_mat_t)A.data[jj - 1];
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

template <typename idx_t, typename coeff_t>
void apply(CSRMatrix<idx_t, coeff_t> const &spmat,
           arma::Mat<coeff_t> const &mat_in, arma::Mat<coeff_t> &mat_out) try {
  apply<idx_t, coeff_t, coeff_t>(spmat, mat_in, mat_out);
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

template <typename idx_t>
void apply(CSRMatrix<idx_t, double> const &spmat,
           arma::Mat<complex> const &mat_in, arma::Mat<complex> &mat_out) try {
  apply<idx_t, double, complex>(spmat, mat_in, mat_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply(CSRMatrix<int32_t, double> const &, arma::cx_mat const &,
                    arma::cx_mat &);
template void apply(CSRMatrix<int64_t, double> const &, arma::cx_mat const &,
                    arma::cx_mat &);

} // namespace xdiag
