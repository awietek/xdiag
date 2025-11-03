// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csr_matrix.hpp"

#include <xdiag/algebra/sparse/csr_matrix_generate.hpp>
#include <xdiag/operators/logic/block.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag {

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

template <typename idx_t, typename coeff_t>
CSRMatrix<idx_t, coeff_t> csr_matrix(OpSum const &op, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  using namespace algebra;
  return std::visit(
      overload{
          [&](Spinhalf const &b1, Spinhalf const &b2) {
            return csr_matrix_generate<idx_t, coeff_t>(op, b1, b2, i0, false);
          },
          [&](tJ const &b1, tJ const &b2) {
            return csr_matrix_generate<idx_t, coeff_t>(op, b1, b2, i0, false);
          },
          [&](Electron const &b1, Electron const &b2) {
            return csr_matrix_generate<idx_t, coeff_t>(op, b1, b2, i0, false);
          },
#ifdef XDIAG_USE_MPI
          [&](SpinhalfDistributed const &, SpinhalfDistributed const &) {
            XDIAG_THROW("Sparse CSRMatrix creation not implemented for "
                        "SpinhalfDistributed blocks");
            return CSRMatrix<idx_t, coeff_t>();
          },
          [&](tJDistributed const &, tJDistributed const &) {
            XDIAG_THROW("Sparse CSRMatrix creation not implemented for "
                        "tJDistributed blocks");
            return CSRMatrix<idx_t, coeff_t>();
          },
          [&](ElectronDistributed const &, ElectronDistributed const &) {
            XDIAG_THROW("Sparse CSRMatrix creation not implemented for "
                        "ElectronDistributed blocks");
            return CSRMatrix<idx_t, coeff_t>();
          },
#endif
          [&](auto &&, auto &&) {
            XDIAG_THROW(fmt::format("Invalid combination of Block types"));
            return CSRMatrix<idx_t, coeff_t>();
          }},
      block_in, block_out);
}
XDIAG_CATCH

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
      idx_t start = csr_mat.rowptr[row];
      idx_t end = csr_mat.rowptr[row + 1];
      for (idx_t i = start; i < end; ++i) {
        idx_t col = csr_mat.col[i];
        coeff_t data = csr_mat.data[i];
        m(row, col) += data;
      }
    }
  } else if (csr_mat.i0 == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (idx_t row = 0; row < csr_mat.nrows; ++row) {
      idx_t start = csr_mat.rowptr[row] - 1;
      idx_t end = csr_mat.rowptr[row + 1] - 1;
      for (idx_t i = start; i < end; ++i) {
        idx_t col = csr_mat.col[i] - 1;
        coeff_t data = csr_mat.data[i];
        m(row, col) += data;
      }
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
