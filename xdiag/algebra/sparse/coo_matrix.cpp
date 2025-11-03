// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix.hpp"

#include <xdiag/algebra/sparse/coo_matrix_generate.hpp>
#include <xdiag/operators/logic/block.hpp>

namespace xdiag {

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

template <typename idx_t, typename coeff_t>
COOMatrix<idx_t, coeff_t> coo_matrix(OpSum const &op, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  using namespace algebra;

  return std::visit(
      overload{[&](Spinhalf const &b1, Spinhalf const &b2) {
                 return coo_matrix_generate<idx_t, coeff_t>(op, b1, b2, i0);
               },
               [&](tJ const &b1, tJ const &b2) {
                 return coo_matrix_generate<idx_t, coeff_t>(op, b1, b2, i0);
               },
               [&](Electron const &b1, Electron const &b2) {
                 return coo_matrix_generate<idx_t, coeff_t>(op, b1, b2, i0);
               },
#ifdef XDIAG_USE_MPI
               [&](SpinhalfDistributed const &, SpinhalfDistributed const &) {
                 XDIAG_THROW("Sparse COOMatrix creation not implemented for "
                             "SpinhalfDistributed blocks");
                 return COOMatrix<idx_t, coeff_t>();
               },
               [&](tJDistributed const &, tJDistributed const &) {
                 XDIAG_THROW("Sparse COOMatrix creation not implemented for "
                             "tJDistributed blocks");
                 return COOMatrix<idx_t, coeff_t>();
               },
               [&](ElectronDistributed const &, ElectronDistributed const &) {
                 XDIAG_THROW("Sparse COOMatrix creation not implemented for "
                             "ElectronDistributed blocks");
                 return COOMatrix<idx_t, coeff_t>();
               },
#endif
               [&](auto &&, auto &&) {
                 XDIAG_THROW(fmt::format("Invalid combination of Block types"));
                 return COOMatrix<idx_t, coeff_t>();
               }},
      block_in, block_out);
}
XDIAG_CATCH

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
