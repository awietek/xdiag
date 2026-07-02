// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csr_matrix.hpp"

#include <algorithm>
#include <numeric>

#include <xdiag/algebra/ishermitian.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/kernels/sparse/sparse_build.hpp>
#include <xdiag/kernels/sparse/valid.hpp>
#include <xdiag/operators/hc.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

// Dispatch overview
// -----------------
// csr_matrix routes through the same two-layer dispatch as matrix():
//
// Layer 1 — Block type (variant dispatch via visit_same_type):
//   csr_matrix(OpSum, Block, Block, idx_t)
//     visit_same_type unwraps both Block variants to their concrete type
//     and calls csr_matrix_impl(OpSum, ConcreteBlock, ConcreteBlock, idx_t).
//
// Layer 2 — csr_matrix_impl<block_t>: one body for every block type. The
//     block's dispatch_basis overload (kernels/blocks/<block>/dispatch_basis
//     .hpp) resolves the shared_ptr<Basis> to a concrete BasisOnTheFly<...>
//     via a type-id table. The numerical kernel is selected by block_t in
//     kernels.cpp. A block type with no dispatch_basis overload is a compile
//     error here, never a silent fallback. The basis-level two-pass build lives
//     in csr_build:
//       Pass 1: csr_matrix_nnz  — counts nonzeros per row
//       Alloc:  builds rowptr and offset from per-row counts; allocates
//               col and data arrays
//       Pass 2: csr_matrix_fill — fills col/data via atomic slot assignment
//       Post:   sorts each row's entries by column index
//
// Kernels — kernels::csr_matrix_nnz / csr_matrix_fill<block_t> (kernels.cpp),
//     instantiated per (block_t, basis type). Offset computation, rowptr
//     construction, and sorting are done here (block level) so they need not be
//     repeated for new block types.

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag {

// Layer 2: block-generic orchestration. The per-basis-type build (dispatch +
// two-pass fill + sort) lives in the shared build_csr_arrays (sparse_build.cpp),
// called here with transpose=false for CSR (groups = rows, ndim = nrows).
template <typename idx_t, typename coeff_t, typename block_t>
static CSRMatrix<idx_t, coeff_t>
csr_matrix_impl(OpSum const &ops, block_t const &block_in,
                block_t const &block_out, idx_t i0) try {
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);
  arma::Col<idx_t> rowptr, col;
  arma::Col<coeff_t> data;
  build_csr_arrays<idx_t, coeff_t, block_t>(ops, block_in, block_out, nrows, i0,
                                            /*transpose=*/false, rowptr, col,
                                            data);
  bool isherm = ishermitian(ops, block_in);
  return CSRMatrix<idx_t, coeff_t>{nrows, ncols, rowptr, col, data, i0, isherm};
}
XDIAG_CATCH

// Layer 1: unwrap Block variant, then call the block-generic Layer 2.
template <typename idx_t, typename coeff_t>
CSRMatrix<idx_t, coeff_t> csr_matrix(OpSum const &ops, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  CSRMatrix<idx_t, coeff_t> result;
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
          result = csr_matrix_impl<idx_t, coeff_t>(ops, bin, bout, i0);
        }
      },
      "Type mismatch of Block types");
  return result;
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
