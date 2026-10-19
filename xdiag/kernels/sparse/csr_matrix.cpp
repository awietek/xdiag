// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csr_matrix.hpp"

#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/compile.hpp>
#include <xdiag/kernels/sparse/sparse_build.hpp>
#include <xdiag/kernels/sparse/valid.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag {

// csr_matrix compiles ops (see kernels/compile.hpp) and builds the arrays in
// the two passes of sparse_build.cpp: per-row counts, then fill and sort.
template <typename idx_t, typename coeff_t>
CSRMatrix<idx_t, coeff_t> csr_matrix(OpSum const &ops, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  check_valid_sparse_matrix<idx_t, coeff_t>(ops, block_in, block_out, i0);
  auto ops_compiled = kernels::compile(ops, block_in, block_out);
  auto counts = sparse_counts<coeff_t>(ops_compiled, /*transpose=*/false);
  int64_t nnz = std::accumulate(counts.begin(), counts.end(), (int64_t)0);
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);
  arma::Col<idx_t> rowptr(nrows + 1), col(nnz);
  arma::Col<coeff_t> data(nnz);
  sparse_fill(ops_compiled, counts, i0, /*transpose=*/false, rowptr.memptr(),
              col.memptr(), data.memptr());
  bool isherm = ishermitian(ops, block_in);
  return CSRMatrix<idx_t, coeff_t>{nrows, ncols, rowptr, col, data, i0, isherm};
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
CSRMatrix<idx_t, coeff_t> csr_matrix(OpSum const &ops, Block const &blocki,
                                     idx_t i0) try {
  auto blocko = block(ops, blocki);
  return csr_matrix<idx_t, coeff_t>(ops, blocki, blocko, i0);
}
XDIAG_CATCH

template CSRMatrix<int32_t, double>
csr_matrix<int32_t, double>(OpSum const &, Block const &, int32_t);
template CSRMatrix<int32_t, complex>
csr_matrix<int32_t, complex>(OpSum const &, Block const &, int32_t);
template CSRMatrix<int64_t, double>
csr_matrix<int64_t, double>(OpSum const &, Block const &, int64_t);
template CSRMatrix<int64_t, complex>
csr_matrix<int64_t, complex>(OpSum const &, Block const &, int64_t);

template CSRMatrix<int32_t, double> csr_matrix<int32_t, double>(OpSum const &,
                                                                Block const &,
                                                                Block const &,
                                                                int32_t);
template CSRMatrix<int32_t, complex> csr_matrix<int32_t, complex>(OpSum const &,
                                                                  Block const &,
                                                                  Block const &,
                                                                  int32_t);
template CSRMatrix<int64_t, double> csr_matrix<int64_t, double>(OpSum const &,
                                                                Block const &,
                                                                Block const &,
                                                                int64_t);
template CSRMatrix<int64_t, complex> csr_matrix<int64_t, complex>(OpSum const &,
                                                                  Block const &,
                                                                  Block const &,
                                                                  int64_t);

// Named convenience wrappers.
CSRMatrix<int64_t, double> csr_matrix(OpSum const &ops, Block const &block,
                                      int64_t i0) try {
  return csr_matrix<int64_t, double>(ops, block, i0);
}
XDIAG_CATCH

CSRMatrix<int64_t, double> csr_matrix(OpSum const &ops, Block const &block_in,
                                      Block const &block_out, int64_t i0) try {
  return csr_matrix<int64_t, double>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

CSRMatrix<int64_t, complex> csr_matrixC(OpSum const &ops, Block const &block,
                                        int64_t i0) try {
  return csr_matrix<int64_t, complex>(ops, block, i0);
}
XDIAG_CATCH

CSRMatrix<int64_t, complex> csr_matrixC(OpSum const &ops, Block const &block_in,
                                        Block const &block_out,
                                        int64_t i0) try {
  return csr_matrix<int64_t, complex>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

CSRMatrix<int32_t, double> csr_matrix_32(OpSum const &ops, Block const &block,
                                         int32_t i0) try {
  return csr_matrix<int32_t, double>(ops, block, i0);
}
XDIAG_CATCH

CSRMatrix<int32_t, double> csr_matrix_32(OpSum const &ops,
                                         Block const &block_in,
                                         Block const &block_out,
                                         int32_t i0) try {
  return csr_matrix<int32_t, double>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

CSRMatrix<int32_t, complex> csr_matrixC_32(OpSum const &ops, Block const &block,
                                           int32_t i0) try {
  return csr_matrix<int32_t, complex>(ops, block, i0);
}
XDIAG_CATCH

CSRMatrix<int32_t, complex> csr_matrixC_32(OpSum const &ops,
                                           Block const &block_in,
                                           Block const &block_out,
                                           int32_t i0) try {
  return csr_matrix<int32_t, complex>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

// Two-phase build into caller-owned storage (Julia wrapper). Phase 1: the
// per-row (per-column with transpose) nonzero counts.
template <typename coeff_t>
std::vector<int64_t> csr_matrix_nnz(OpSum const &ops, Block const &block_in,
                                    Block const &block_out,
                                    bool transpose) try {
  return sparse_counts<coeff_t>(kernels::compile(ops, block_in, block_out),
                                transpose);
}
XDIAG_CATCH

// Phase 2: fill the caller's rowptr/col/data, using the phase-1 counts.
template <typename idx_t, typename coeff_t>
void csr_matrix_fill(OpSum const &ops, Block const &block_in,
                     Block const &block_out,
                     std::vector<int64_t> const &n_elements_in_row,
                     idx_t *rowptr, idx_t *col, coeff_t *data, idx_t i0,
                     bool transpose) try {
  check_valid_sparse_matrix<idx_t, coeff_t>(ops, block_in, block_out, i0);
  sparse_fill(kernels::compile(ops, block_in, block_out), n_elements_in_row, i0,
              transpose, rowptr, col, data);
}
XDIAG_CATCH

template std::vector<int64_t>
csr_matrix_nnz<double>(OpSum const &, Block const &, Block const &, bool);
template std::vector<int64_t>
csr_matrix_nnz<complex>(OpSum const &, Block const &, Block const &, bool);

template void csr_matrix_fill<int32_t, double>(OpSum const &, Block const &,
                                               Block const &,
                                               std::vector<int64_t> const &,
                                               int32_t *, int32_t *, double *,
                                               int32_t, bool);
template void csr_matrix_fill<int64_t, double>(OpSum const &, Block const &,
                                               Block const &,
                                               std::vector<int64_t> const &,
                                               int64_t *, int64_t *, double *,
                                               int64_t, bool);
template void csr_matrix_fill<int32_t, complex>(OpSum const &, Block const &,
                                                Block const &,
                                                std::vector<int64_t> const &,
                                                int32_t *, int32_t *, complex *,
                                                int32_t, bool);
template void csr_matrix_fill<int64_t, complex>(OpSum const &, Block const &,
                                                Block const &,
                                                std::vector<int64_t> const &,
                                                int64_t *, int64_t *, complex *,
                                                int64_t, bool);

template <typename idx_t, typename coeff_t>
arma::Mat<coeff_t> to_dense(CSRMatrix<idx_t, coeff_t> const &csr_mat) try {
  int64_t nnz = csr_mat.data.size();
  if ((csr_mat.nrows == 0) || (csr_mat.ncols == 0)) {
    return arma::Mat<coeff_t>();
  }

  if (csr_mat.rowptr.size() != csr_mat.nrows + 1) {
    XDIAG_THROW(fmt::format(
        "Number of rowptr entries ({}) does not match number of rows (+1) ({})",
        csr_mat.rowptr.size(), csr_mat.nrows + 1));
  }
  if (csr_mat.col.size() != nnz) {
    XDIAG_THROW(
        "Number of column entries does not match number of data entries");
  }

  arma::Mat<coeff_t> m(csr_mat.nrows, csr_mat.ncols, arma::fill::zeros);

  if (csr_mat.i0 == 0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (idx_t row = 0; row < csr_mat.nrows; ++row) {
      for (idx_t i = csr_mat.rowptr[row]; i < csr_mat.rowptr[row + 1]; ++i)
        m(row, csr_mat.col[i]) += csr_mat.data[i];
    }
  } else if (csr_mat.i0 == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (idx_t row = 0; row < csr_mat.nrows; ++row) {
      for (idx_t i = csr_mat.rowptr[row] - 1; i < csr_mat.rowptr[row + 1] - 1;
           ++i)
        m(row, csr_mat.col[i] - 1) += csr_mat.data[i];
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Invalid zero index i0. Must be either 0 or 1, but got i0={}",
        csr_mat.i0));
  }
  return m;
}
XDIAG_CATCH

template arma::mat to_dense(CSRMatrix<int32_t, double> const &);
template arma::mat to_dense(CSRMatrix<int64_t, double> const &);
template arma::cx_mat to_dense(CSRMatrix<int32_t, complex> const &);
template arma::cx_mat to_dense(CSRMatrix<int64_t, complex> const &);

} // namespace xdiag
