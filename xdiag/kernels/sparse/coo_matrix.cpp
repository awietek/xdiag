// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix.hpp"

#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/compile.hpp>
#include <xdiag/kernels/elements.hpp>
#include <xdiag/kernels/sparse/sparse_build.hpp>
#include <xdiag/kernels/sparse/valid.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag {

// COO entries grouped by row, like CSR: the per-row counts of sparse_build.cpp
// give every row its slots, and each row is filled by one thread in a fixed
// order. The result is row-major and does not depend on the thread count.
template <typename idx_t, typename coeff_t>
static void coo_fill(kernels::CompiledOpSum const &ops,
                     std::vector<int64_t> const &counts, idx_t i0, idx_t *row,
                     idx_t *col, coeff_t *data) {
  std::vector<int64_t> offset(counts.size() + 1, 0);
  std::partial_sum(counts.begin(), counts.end(), offset.begin() + 1);
  kernels::elements<coeff_t>(ops, [&](int64_t r, int64_t c, coeff_t val) {
    int64_t k = offset[r]++;
    row[k] = (idx_t)(r + i0);
    col[k] = (idx_t)(c + i0);
    data[k] = val;
  });
}

template <typename idx_t, typename coeff_t>
COOMatrix<idx_t, coeff_t> coo_matrix(OpSum const &ops, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  check_valid_sparse_matrix<idx_t, coeff_t>(ops, block_in, block_out, i0);
  auto ops_compiled = kernels::compile(ops, block_in, block_out);
  auto counts = sparse_counts<coeff_t>(ops_compiled, /*transpose=*/false);
  int64_t nnz = std::accumulate(counts.begin(), counts.end(), (int64_t)0);
  arma::Col<idx_t> row(nnz), col(nnz);
  arma::Col<coeff_t> data(nnz);
  coo_fill(ops_compiled, counts, i0, row.memptr(), col.memptr(), data.memptr());
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);
  bool isherm = ishermitian(ops, block_in);
  return COOMatrix<idx_t, coeff_t>{nrows, ncols, row, col, data, i0, isherm};
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
COOMatrix<idx_t, coeff_t> coo_matrix(OpSum const &ops, Block const &blocki,
                                     idx_t i0) try {
  auto blocko = block(ops, blocki);
  return coo_matrix<idx_t, coeff_t>(ops, blocki, blocko, i0);
}
XDIAG_CATCH

template COOMatrix<int32_t, double>
coo_matrix<int32_t, double>(OpSum const &, Block const &, int32_t);
template COOMatrix<int32_t, complex>
coo_matrix<int32_t, complex>(OpSum const &, Block const &, int32_t);
template COOMatrix<int64_t, double>
coo_matrix<int64_t, double>(OpSum const &, Block const &, int64_t);
template COOMatrix<int64_t, complex>
coo_matrix<int64_t, complex>(OpSum const &, Block const &, int64_t);

template COOMatrix<int32_t, double> coo_matrix<int32_t, double>(OpSum const &,
                                                                Block const &,
                                                                Block const &,
                                                                int32_t);
template COOMatrix<int32_t, complex> coo_matrix<int32_t, complex>(OpSum const &,
                                                                  Block const &,
                                                                  Block const &,
                                                                  int32_t);
template COOMatrix<int64_t, double> coo_matrix<int64_t, double>(OpSum const &,
                                                                Block const &,
                                                                Block const &,
                                                                int64_t);
template COOMatrix<int64_t, complex> coo_matrix<int64_t, complex>(OpSum const &,
                                                                  Block const &,
                                                                  Block const &,
                                                                  int64_t);

// Named convenience wrappers (no template syntax at call sites).
COOMatrix<int64_t, double> coo_matrix(OpSum const &ops, Block const &block,
                                      int64_t i0) try {
  return coo_matrix<int64_t, double>(ops, block, i0);
}
XDIAG_CATCH

COOMatrix<int64_t, double> coo_matrix(OpSum const &ops, Block const &block_in,
                                      Block const &block_out, int64_t i0) try {
  return coo_matrix<int64_t, double>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

COOMatrix<int64_t, complex> coo_matrixC(OpSum const &ops, Block const &block,
                                        int64_t i0) try {
  return coo_matrix<int64_t, complex>(ops, block, i0);
}
XDIAG_CATCH

COOMatrix<int64_t, complex> coo_matrixC(OpSum const &ops, Block const &block_in,
                                        Block const &block_out,
                                        int64_t i0) try {
  return coo_matrix<int64_t, complex>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

COOMatrix<int32_t, double> coo_matrix_32(OpSum const &ops, Block const &block,
                                         int32_t i0) try {
  return coo_matrix<int32_t, double>(ops, block, i0);
}
XDIAG_CATCH

COOMatrix<int32_t, double> coo_matrix_32(OpSum const &ops,
                                         Block const &block_in,
                                         Block const &block_out,
                                         int32_t i0) try {
  return coo_matrix<int32_t, double>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

COOMatrix<int32_t, complex> coo_matrixC_32(OpSum const &ops, Block const &block,
                                           int32_t i0) try {
  return coo_matrix<int32_t, complex>(ops, block, i0);
}
XDIAG_CATCH

COOMatrix<int32_t, complex> coo_matrixC_32(OpSum const &ops,
                                           Block const &block_in,
                                           Block const &block_out,
                                           int32_t i0) try {
  return coo_matrix<int32_t, complex>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

// Two-phase build into caller-owned storage (Julia wrapper). Phase 1 returns
// the total nnz, phase 2 fills the caller's row/col/data.
template <typename coeff_t>
int64_t coo_matrix_nnz(OpSum const &ops, Block const &block_in,
                       Block const &block_out) try {
  auto counts = sparse_counts<coeff_t>(
      kernels::compile(ops, block_in, block_out), /*transpose=*/false);
  return std::accumulate(counts.begin(), counts.end(), (int64_t)0);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
void coo_matrix_fill(OpSum const &ops, Block const &block_in,
                     Block const &block_out, int64_t nnz_capacity, idx_t *row,
                     idx_t *col, coeff_t *data, idx_t i0) try {
  check_valid_sparse_matrix<idx_t, coeff_t>(ops, block_in, block_out, i0);
  auto ops_compiled = kernels::compile(ops, block_in, block_out);
  auto counts = sparse_counts<coeff_t>(ops_compiled, /*transpose=*/false);
  int64_t nnz = std::accumulate(counts.begin(), counts.end(), (int64_t)0);
  if (nnz > nnz_capacity) {
    XDIAG_THROW("coo_matrix_fill: the number of nonzeros exceeds the "
                "allocated capacity");
  }
  coo_fill(ops_compiled, counts, i0, row, col, data);
}
XDIAG_CATCH

template int64_t coo_matrix_nnz<double>(OpSum const &, Block const &,
                                        Block const &);
template int64_t coo_matrix_nnz<complex>(OpSum const &, Block const &,
                                         Block const &);
template void coo_matrix_fill<int32_t, double>(OpSum const &, Block const &,
                                               Block const &, int64_t,
                                               int32_t *, int32_t *, double *,
                                               int32_t);
template void coo_matrix_fill<int64_t, double>(OpSum const &, Block const &,
                                               Block const &, int64_t,
                                               int64_t *, int64_t *, double *,
                                               int64_t);
template void coo_matrix_fill<int32_t, complex>(OpSum const &, Block const &,
                                                Block const &, int64_t,
                                                int32_t *, int32_t *, complex *,
                                                int32_t);
template void coo_matrix_fill<int64_t, complex>(OpSum const &, Block const &,
                                                Block const &, int64_t,
                                                int64_t *, int64_t *, complex *,
                                                int64_t);

template <typename idx_t, typename coeff_t>
arma::Mat<coeff_t> to_dense(COOMatrix<idx_t, coeff_t> const &coo_mat) try {
  if ((coo_mat.nrows == 0) || (coo_mat.ncols == 0)) {
    return arma::Mat<coeff_t>();
  }

  int64_t nnz = coo_mat.data.size();
  if (coo_mat.row.size() != nnz) {
    XDIAG_THROW("Number of row entries does not match number of data entries");
  }
  if (coo_mat.col.size() != nnz) {
    XDIAG_THROW(
        "Number of column entries does not match number of data entries");
  }

  arma::Mat<coeff_t> m(coo_mat.nrows, coo_mat.ncols, arma::fill::zeros);

  if (coo_mat.i0 == 0) {
    for (int64_t idx = 0; idx < nnz; ++idx) {
      int64_t i = coo_mat.row[idx];
      int64_t j = coo_mat.col[idx];
      coeff_t val = coo_mat.data[idx];
      m(i, j) += val;
    }
  } else if (coo_mat.i0 == 1) {
    for (int64_t idx = 0; idx < nnz; ++idx) {
      int64_t i = coo_mat.row[idx] - 1;
      int64_t j = coo_mat.col[idx] - 1;
      coeff_t val = coo_mat.data[idx];
      m(i, j) += val;
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Invalid zero index i0. Must be either 0 or 1, but got i0={}",
        coo_mat.i0));
  }
  return m;
}
XDIAG_CATCH

template arma::mat to_dense(COOMatrix<int32_t, double> const &);
template arma::mat to_dense(COOMatrix<int64_t, double> const &);
template arma::cx_mat to_dense(COOMatrix<int32_t, complex> const &);
template arma::cx_mat to_dense(COOMatrix<int64_t, complex> const &);

} // namespace xdiag
