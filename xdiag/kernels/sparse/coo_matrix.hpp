// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/sparse/sparse_matrix_types.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

// int64_t, double
XDIAG_API COOMatrix<int64_t, double>
coo_matrix(OpSum const &ops, Block const &block, int64_t i0 = 0);
XDIAG_API COOMatrix<int64_t, double> coo_matrix(OpSum const &ops,
                                                Block const &block_in,
                                                Block const &block_out,
                                                int64_t i0 = 0);
// int64_t, complex
XDIAG_API COOMatrix<int64_t, complex>
coo_matrixC(OpSum const &ops, Block const &block, int64_t i0 = 0);
XDIAG_API COOMatrix<int64_t, complex> coo_matrixC(OpSum const &ops,
                                                  Block const &block_in,
                                                  Block const &block_out,
                                                  int64_t i0 = 0);
// int32_t, double
XDIAG_API COOMatrix<int32_t, double>
coo_matrix_32(OpSum const &ops, Block const &block, int32_t i0 = 0);
XDIAG_API COOMatrix<int32_t, double> coo_matrix_32(OpSum const &ops,
                                                   Block const &block_in,
                                                   Block const &block_out,
                                                   int32_t i0 = 0);
// int32_t, complex
XDIAG_API COOMatrix<int32_t, complex>
coo_matrixC_32(OpSum const &ops, Block const &block, int32_t i0 = 0);
XDIAG_API COOMatrix<int32_t, complex> coo_matrixC_32(OpSum const &ops,
                                                     Block const &block_in,
                                                     Block const &block_out,
                                                     int32_t i0 = 0);
template <typename idx_t, typename coeff_t>
XDIAG_API COOMatrix<idx_t, coeff_t>
coo_matrix(OpSum const &ops, Block const &block, idx_t i0 = 0);

template <typename idx_t, typename coeff_t>
XDIAG_API COOMatrix<idx_t, coeff_t>
coo_matrix(OpSum const &ops, Block const &block_in, Block const &block_out,
           idx_t i0 = 0);

template <typename idx_t, typename coeff_t>
XDIAG_API arma::Mat<coeff_t> to_dense(COOMatrix<idx_t, coeff_t> const &coo_mat);

// Two-phase build for caller-allocated storage (used by the Julia wrapper,
// which owns the row/col/data arrays).
//
// Phase 1 — coo_matrix_nnz: returns the total number of nonzeros. The caller
// allocates row (nnz), col (nnz) and data (nnz).
//
// Phase 2 — coo_matrix_fill: populates the caller's row/col/data. nnz_capacity
// is the length the caller allocated (the phase-1 result); the routine throws
// rather than overflow if the recomputed count exceeds it (a guard against the
// OpenMP thread count changing between the two calls). i0 is the index base
// (0 or 1). No pointer is retained past the call. coeff_t in {double, complex};
// idx_t in {int32_t, int64_t}. Distributed blocks throw.
template <typename coeff_t>
XDIAG_API int64_t coo_matrix_nnz(OpSum const &ops, Block const &block_in,
                                 Block const &block_out);

template <typename idx_t, typename coeff_t>
XDIAG_API void coo_matrix_fill(OpSum const &ops, Block const &block_in,
                               Block const &block_out, int64_t nnz_capacity,
                               idx_t *row, idx_t *col, coeff_t *data, idx_t i0);
} // namespace xdiag
