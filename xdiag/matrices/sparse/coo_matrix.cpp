// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix.hpp"

#include <numeric>

#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/matrices/kernels.hpp>
#include <xdiag/matrices/sparse/valid.hpp>
#include <xdiag/matrices/spinhalf/dispatch_basis.hpp>
#include <xdiag/matrices/spinhalf/matrix_policy.hpp>
#include <xdiag/operators/hc.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

// Dispatch overview
// -----------------
// coo_matrix routes through the same two-layer dispatch as matrix():
//
// Layer 1 — Block type (variant dispatch via visit_same_type):
//   coo_matrix(OpSum, Block, Block, idx_t)
//     visit_same_type unwraps both Block variants to their concrete type
//     and calls coo_matrix(OpSum, ConcreteBlock, ConcreteBlock, idx_t).
//
// Layer 2 — Basis type (runtime dispatch via dispatch_basis):
//   coo_matrix(OpSum, Spinhalf, Spinhalf, idx_t)
//     dispatch_basis resolves the shared_ptr<Basis> to a concrete
//     BasisOnTheFly<...> type, then orchestrates the two-pass COO build:
//       Pass 1: coo_matrix_nnz  — counts NNZ (per thread with OMP)
//       Alloc:  allocates row/col/data arrays from the NNZ count
//       Pass 2: coo_matrix_fill — fills the pre-allocated arrays
//
// Kernels — spinhalf::coo_matrix_nnz / coo_matrix_fill (kernels.cpp)
//     Defined in spinhalf/kernels.cpp; instantiated per basis type via
//     the instantiation-group mechanism.  The two-pass split keeps the
//     memory allocation at the block level so it is not repeated when
//     new block types are added.

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag {

// Layer 2: Spinhalf-specific function.  Orchestrates the two-pass COO build
// using dispatch_basis to obtain concrete BasisOnTheFly objects.
template <typename idx_t, typename coeff_t>
static COOMatrix<idx_t, coeff_t>
coo_matrix(OpSum const &ops, Spinhalf const &block_in,
           Spinhalf const &block_out, idx_t i0) try {
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);
  arma::Col<idx_t> rows, cols;
  arma::Col<coeff_t> data;

  matrices::spinhalf::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
#ifdef _OPENMP
        // Pass 1: count NNZ per thread.
        auto nnz_thread =
            matrices::coo_matrix_nnz<matrices::spinhalf::MatrixPolicy, coeff_t>(
                ops, basis_in, basis_out);
        int64_t nnz =
            std::accumulate(nnz_thread.begin(), nnz_thread.end(), (int64_t)0);
        rows.resize(nnz);
        cols.resize(nnz);
        data.resize(nnz);
        // Pass 2: fill using per-thread offsets.
        matrices::coo_matrix_fill<matrices::spinhalf::MatrixPolicy, coeff_t>(
            ops, basis_in, basis_out, nnz_thread, rows.memptr(), cols.memptr(),
            data.memptr(), i0);
#else
        // Pass 1: count total NNZ.
        int64_t nnz =
            matrices::coo_matrix_nnz<matrices::spinhalf::MatrixPolicy, coeff_t>(
                ops, basis_in, basis_out);
        rows.resize(nnz);
        cols.resize(nnz);
        data.resize(nnz);
        // Pass 2: fill sequentially.
        matrices::coo_matrix_fill<matrices::spinhalf::MatrixPolicy, coeff_t>(
            ops, basis_in, basis_out, rows.memptr(), cols.memptr(),
            data.memptr(), i0);
#endif
      });

  bool isherm = ishermitian(ops, block_in);
  return COOMatrix<idx_t, coeff_t>{nrows, ncols, rows, cols, data, i0, isherm};
}
XDIAG_CATCH

// Layer 1: unwrap Block variant; validate with concrete types, then call the
// block-specific Layer 2 function (currently Spinhalf only).
template <typename idx_t, typename coeff_t>
COOMatrix<idx_t, coeff_t> coo_matrix(OpSum const &ops, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  COOMatrix<idx_t, coeff_t> result;
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        check_valid_sparse_matrix<idx_t, coeff_t>(ops, bin, bout, i0);
        result = coo_matrix<idx_t, coeff_t>(ops, bin, bout, i0);
      },
      "Type mismatch of Block types");
  return result;
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
