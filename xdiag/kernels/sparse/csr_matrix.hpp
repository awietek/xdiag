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
XDIAG_API CSRMatrix<int64_t, double>
csr_matrix(OpSum const &ops, Block const &block, int64_t i0 = 0);
XDIAG_API CSRMatrix<int64_t, double> csr_matrix(OpSum const &ops,
                                                Block const &block_in,
                                                Block const &block_out,
                                                int64_t i0 = 0);

// int64_t, complex
XDIAG_API CSRMatrix<int64_t, complex>
csr_matrixC(OpSum const &ops, Block const &block, int64_t i0 = 0);
XDIAG_API CSRMatrix<int64_t, complex> csr_matrixC(OpSum const &ops,
                                                  Block const &block_in,
                                                  Block const &block_out,
                                                  int64_t i0 = 0);
// int32_t, double
XDIAG_API CSRMatrix<int32_t, double>
csr_matrix_32(OpSum const &ops, Block const &block, int32_t i0 = 0);
XDIAG_API CSRMatrix<int32_t, double> csr_matrix_32(OpSum const &ops,
                                                   Block const &block_in,
                                                   Block const &block_out,
                                                   int32_t i0 = 0);
// int32_t, complex
XDIAG_API CSRMatrix<int32_t, complex>
csr_matrixC_32(OpSum const &ops, Block const &block, int32_t i0 = 0);
XDIAG_API CSRMatrix<int32_t, complex> csr_matrixC_32(OpSum const &ops,
                                                     Block const &block_in,
                                                     Block const &block_out,
                                                     int32_t i0 = 0);

template <typename idx_t, typename coeff_t>
XDIAG_API CSRMatrix<idx_t, coeff_t>
csr_matrix(OpSum const &ops, Block const &block, idx_t i0 = 0);

template <typename idx_t, typename coeff_t>
XDIAG_API CSRMatrix<idx_t, coeff_t>
csr_matrix(OpSum const &ops, Block const &block_in, Block const &block_out,
           idx_t i0 = 0);

template <typename idx_t, typename coeff_t>
XDIAG_API arma::Mat<coeff_t> to_dense(CSRMatrix<idx_t, coeff_t> const &csr_mat);

// Two-phase build for caller-allocated storage (used by the Julia wrapper,
// which owns the rowptr/col/data arrays).
//
// Phase 1 — csr_matrix_nnz: returns the per-row nonzero counts (length
// size(block_out)). The caller sums these to obtain nnz, then allocates
// rowptr (size(block_out)+1), col (nnz) and data (nnz).
//
// Phase 2 — csr_matrix_fill: populates the caller's rowptr/col/data. The
// phase-1 counts must be passed back in as n_elements_in_row (used to build
// rowptr and the per-row write offsets, avoiding a second counting pass). i0
// is the index base (0 or 1). Entries within each row are sorted by column.
// No pointer is retained past the call; the buffers must outlive it only.
// coeff_t in {double, complex}; idx_t in {int32_t, int64_t}. Distributed
// blocks throw.
template <typename coeff_t>
XDIAG_API std::vector<int64_t>
csr_matrix_nnz(OpSum const &ops, Block const &block_in, Block const &block_out);

template <typename idx_t, typename coeff_t>
XDIAG_API void csr_matrix_fill(OpSum const &ops, Block const &block_in,
                               Block const &block_out,
                               std::vector<int64_t> const &n_elements_in_row,
                               idx_t *rowptr, idx_t *col, coeff_t *data,
                               idx_t i0);
} // namespace xdiag
