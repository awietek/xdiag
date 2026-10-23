// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csc_matrix.hpp"

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

// csc_matrix compiles ops (see kernels/compile.hpp) and builds the arrays in
// the two passes of sparse_build.cpp, grouped by column (transpose).
template <typename idx_t, typename coeff_t>
CSCMatrix<idx_t, coeff_t> csc_matrix(OpSum const &ops, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  check_valid_sparse_matrix<idx_t, coeff_t>(ops, block_in, block_out, i0);
  auto ops_compiled = kernels::compile(ops, block_in, block_out);
  auto counts = sparse_counts<coeff_t>(ops_compiled, /*transpose=*/true);
  int64_t nnz = std::accumulate(counts.begin(), counts.end(), (int64_t)0);
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);
  arma::Col<idx_t> colptr(ncols + 1), row(nnz);
  arma::Col<coeff_t> data(nnz);
  sparse_fill(ops_compiled, counts, i0, /*transpose=*/true, colptr.memptr(),
              row.memptr(), data.memptr());
  bool isherm = ishermitian(ops, block_in);
  return CSCMatrix<idx_t, coeff_t>{nrows, ncols, colptr, row, data, i0, isherm};
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
CSCMatrix<idx_t, coeff_t> csc_matrix(OpSum const &ops, Block const &blocki,
                                     idx_t i0) try {
  auto blocko = block(ops, blocki);
  return csc_matrix<idx_t, coeff_t>(ops, blocki, blocko, i0);
}
XDIAG_CATCH

template CSCMatrix<int32_t, double>
csc_matrix<int32_t, double>(OpSum const &, Block const &, int32_t);
template CSCMatrix<int32_t, complex>
csc_matrix<int32_t, complex>(OpSum const &, Block const &, int32_t);
template CSCMatrix<int64_t, double>
csc_matrix<int64_t, double>(OpSum const &, Block const &, int64_t);
template CSCMatrix<int64_t, complex>
csc_matrix<int64_t, complex>(OpSum const &, Block const &, int64_t);

template CSCMatrix<int32_t, double> csc_matrix<int32_t, double>(OpSum const &,
                                                                Block const &,
                                                                Block const &,
                                                                int32_t);
template CSCMatrix<int32_t, complex> csc_matrix<int32_t, complex>(OpSum const &,
                                                                  Block const &,
                                                                  Block const &,
                                                                  int32_t);
template CSCMatrix<int64_t, double> csc_matrix<int64_t, double>(OpSum const &,
                                                                Block const &,
                                                                Block const &,
                                                                int64_t);
template CSCMatrix<int64_t, complex> csc_matrix<int64_t, complex>(OpSum const &,
                                                                  Block const &,
                                                                  Block const &,
                                                                  int64_t);

// Named convenience wrappers.
CSCMatrix<int64_t, double> csc_matrix(OpSum const &ops, Block const &block,
                                      int64_t i0) try {
  return csc_matrix<int64_t, double>(ops, block, i0);
}
XDIAG_CATCH

CSCMatrix<int64_t, double> csc_matrix(OpSum const &ops, Block const &block_in,
                                      Block const &block_out, int64_t i0) try {
  return csc_matrix<int64_t, double>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

CSCMatrix<int64_t, complex> csc_matrixC(OpSum const &ops, Block const &block,
                                        int64_t i0) try {
  return csc_matrix<int64_t, complex>(ops, block, i0);
}
XDIAG_CATCH

CSCMatrix<int64_t, complex> csc_matrixC(OpSum const &ops, Block const &block_in,
                                        Block const &block_out,
                                        int64_t i0) try {
  return csc_matrix<int64_t, complex>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

CSCMatrix<int32_t, double> csc_matrix_32(OpSum const &ops, Block const &block,
                                         int32_t i0) try {
  return csc_matrix<int32_t, double>(ops, block, i0);
}
XDIAG_CATCH

CSCMatrix<int32_t, double> csc_matrix_32(OpSum const &ops,
                                         Block const &block_in,
                                         Block const &block_out,
                                         int32_t i0) try {
  return csc_matrix<int32_t, double>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

CSCMatrix<int32_t, complex> csc_matrixC_32(OpSum const &ops, Block const &block,
                                           int32_t i0) try {
  return csc_matrix<int32_t, complex>(ops, block, i0);
}
XDIAG_CATCH

CSCMatrix<int32_t, complex> csc_matrixC_32(OpSum const &ops,
                                           Block const &block_in,
                                           Block const &block_out,
                                           int32_t i0) try {
  return csc_matrix<int32_t, complex>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
arma::Mat<coeff_t> to_dense(CSCMatrix<idx_t, coeff_t> const &csc_mat) try {
  if ((csc_mat.nrows == 0) || (csc_mat.ncols == 0)) {
    return arma::Mat<coeff_t>();
  }
  int64_t nnz = csc_mat.data.size();
  if (csc_mat.colptr.size() != csc_mat.ncols + 1) {
    XDIAG_THROW("Number of colptr entries does not match number of cols (+1)");
  }
  if (csc_mat.row.size() != nnz) {
    XDIAG_THROW("Number of row entries does not match number of data entries");
  }

  arma::Mat<coeff_t> m(csc_mat.nrows, csc_mat.ncols, arma::fill::zeros);

  if (csc_mat.i0 == 0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (idx_t col = 0; col < csc_mat.ncols; ++col) {
      for (idx_t i = csc_mat.colptr[col]; i < csc_mat.colptr[col + 1]; ++i)
        m(csc_mat.row[i], col) += csc_mat.data[i];
    }
  } else if (csc_mat.i0 == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (idx_t col = 0; col < csc_mat.ncols; ++col) {
      for (idx_t i = csc_mat.colptr[col] - 1; i < csc_mat.colptr[col + 1] - 1;
           ++i)
        m(csc_mat.row[i] - 1, col) += csc_mat.data[i];
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Invalid zero index i0. Must be either 0 or 1, but got i0={}",
        csc_mat.i0));
  }
  return m;
}
XDIAG_CATCH

template arma::mat to_dense(CSCMatrix<int32_t, double> const &);
template arma::mat to_dense(CSCMatrix<int64_t, double> const &);
template arma::cx_mat to_dense(CSCMatrix<int32_t, complex> const &);
template arma::cx_mat to_dense(CSCMatrix<int64_t, complex> const &);

} // namespace xdiag
