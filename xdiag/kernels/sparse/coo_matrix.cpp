// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix.hpp"

#include <numeric>

#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/kernels/blocks/dispatch_bases.hpp>
#include <xdiag/kernels/kernels.hpp>
#include <xdiag/kernels/sparse/valid.hpp>
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
//     and calls coo_matrix_impl(OpSum, ConcreteBlock, ConcreteBlock, idx_t).
//
// Layer 2 — coo_matrix_impl<block_t>: one body for every block type. The
//     block's dispatch_basis overload (kernels/blocks/<block>/dispatch_basis
//     .hpp) resolves the shared_ptr<Basis> to a concrete BasisOnTheFly<...>
//     via a type-id table. The numerical kernel is selected by block_t in
//     kernels.cpp. A block type with no dispatch_basis overload is a compile
//     error here, never a silent fallback. The basis-level two-pass build lives
//     in coo_build:
//       Pass 1: coo_matrix_nnz  — counts NNZ (per thread with OMP)
//       Alloc:  allocates row/col/data arrays from the NNZ count
//       Pass 2: coo_matrix_fill — fills the pre-allocated arrays
//
// Kernels — kernels::coo_matrix_nnz / coo_matrix_fill<block_t> (kernels.cpp),
//     instantiated per (block_t, basis type) via the instantiation-group
//     mechanism. The two-pass split keeps the memory allocation at the block
//     level so it is not repeated when new block types are added.

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag {

// Basis-level two-pass COO build, generic over the block type and
// concrete BasisOnTheFly type. Shared by the block-generic Layer-2 below.
template <typename idx_t, typename coeff_t, typename block_t, typename basis_t>
static void coo_build(OpSum const &ops, basis_t const &basis_in,
                      basis_t const &basis_out, idx_t i0,
                      arma::Col<idx_t> &rows, arma::Col<idx_t> &cols,
                      arma::Col<coeff_t> &data) {
#ifdef _OPENMP
  // Pass 1: count NNZ per thread.
  auto nnz_thread =
      kernels::coo_matrix_nnz<block_t, coeff_t>(ops, basis_in, basis_out);
  int64_t nnz =
      std::accumulate(nnz_thread.begin(), nnz_thread.end(), (int64_t)0);
  rows.resize(nnz);
  cols.resize(nnz);
  data.resize(nnz);
  // Pass 2: fill using per-thread offsets.
  kernels::coo_matrix_fill<block_t, coeff_t>(ops, basis_in, basis_out,
                                               nnz_thread, rows.memptr(),
                                               cols.memptr(), data.memptr(), i0);
#else
  // Pass 1: count total NNZ.
  int64_t nnz =
      kernels::coo_matrix_nnz<block_t, coeff_t>(ops, basis_in, basis_out);
  rows.resize(nnz);
  cols.resize(nnz);
  data.resize(nnz);
  // Pass 2: fill sequentially.
  kernels::coo_matrix_fill<block_t, coeff_t>(
      ops, basis_in, basis_out, rows.memptr(), cols.memptr(), data.memptr(),
      i0);
#endif
}

// Layer 2: block-generic orchestration. The dispatch_basis overload supplies
// the basis dispatch; a block with no overload is a compile error here, never a
// silent fallback to the Block overload.
template <typename idx_t, typename coeff_t, typename block_t>
static COOMatrix<idx_t, coeff_t>
coo_matrix_impl(OpSum const &ops, block_t const &block_in,
                block_t const &block_out, idx_t i0) try {
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);
  arma::Col<idx_t> rows, cols;
  arma::Col<coeff_t> data;
  kernels::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
        coo_build<idx_t, coeff_t, block_t>(ops, basis_in, basis_out, i0, rows,
                                           cols, data);
      });
  bool isherm = ishermitian(ops, block_in);
  return COOMatrix<idx_t, coeff_t>{nrows, ncols, rows, cols, data, i0, isherm};
}
XDIAG_CATCH

// Layer 1: unwrap Block variant; validate with concrete types, then call the
// block-generic Layer 2.
template <typename idx_t, typename coeff_t>
COOMatrix<idx_t, coeff_t> coo_matrix(OpSum const &ops, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  COOMatrix<idx_t, coeff_t> result;
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        using block_t = std::decay_t<decltype(bin)>;
        if constexpr (is_distributed_v<block_t>) {
          XDIAG_THROW("Cannot build a sparse matrix for a distributed block: "
                      "its Hilbert space is distributed across MPI ranks. Use "
                      "apply(...) instead.");
        } else {
          check_valid_sparse_matrix<idx_t, coeff_t>(ops, bin, bout, i0);
          result = coo_matrix_impl<idx_t, coeff_t>(ops, bin, bout, i0);
        }
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
