// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csc_matrix.hpp"

#include <xdiag/algebra/sparse/csr_matrix_generate.hpp>
#include <xdiag/operators/logic/block.hpp>

namespace xdiag {

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

template <typename idx_t, typename coeff_t>
static CSCMatrix<idx_t, coeff_t> to_csc(CSRMatrix<idx_t, coeff_t> &&csr) {
  return CSCMatrix<idx_t, coeff_t>{
      csr.ncols,          csr.nrows,           std::move(csr.rowptr),
      std::move(csr.col), std::move(csr.data), csr.i0};
}

template <typename idx_t, typename coeff_t>
CSCMatrix<idx_t, coeff_t> csc_matrix(OpSum const &op, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  using namespace algebra;
  return to_csc(std::visit(
      overload{
          [&](Spinhalf const &b1, Spinhalf const &b2) {
            return csr_matrix_generate<idx_t, coeff_t>(op, b1, b2, i0, true);
          },
          [&](tJ const &b1, tJ const &b2) {
            return csr_matrix_generate<idx_t, coeff_t>(op, b1, b2, i0, true);
          },
          [&](Electron const &b1, Electron const &b2) {
            return csr_matrix_generate<idx_t, coeff_t>(op, b1, b2, i0, true);
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
      block_in, block_out));
}
XDIAG_CATCH

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
      idx_t start = csc_mat.colptr[col];
      idx_t end = csc_mat.colptr[col + 1];
      for (idx_t i = start; i < end; ++i) {
        idx_t row = csc_mat.row[i];
        coeff_t data = csc_mat.data[i];
        m(row, col) += data;
      }
    }
  } else if (csc_mat.i0 == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (idx_t col = 0; col < csc_mat.ncols; ++col) {
      idx_t start = csc_mat.colptr[col] - 1;
      idx_t end = csc_mat.colptr[col + 1] - 1;
      for (idx_t i = start; i < end; ++i) {
        idx_t row = csc_mat.row[i] - 1;
        coeff_t data = csc_mat.data[i];
        m(row, col) += data;
      }
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
