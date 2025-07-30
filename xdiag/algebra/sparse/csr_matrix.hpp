// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <xdiag/common.hpp>

#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

XDIAG_API csr_mat csr_matrix(Op const &op, Block const &block, int64_t i0 = 0);
XDIAG_API csr_mat csr_matrix(OpSum const &ops, Block const &block,
                             int64_t i0 = 0);
XDIAG_API csr_cx_mat csr_matrixC(Op const &op, Block const &block,
                                 int64_t i0 = 0);
XDIAG_API csr_cx_mat csr_matrixC(OpSum const &ops, Block const &block,
                                 int64_t i0 = 0);

XDIAG_API csr_mat csr_matrix(Op const &op, Block const &block_in,
                             Block const &block_out, int64_t i0 = 0);
XDIAG_API csr_mat csr_matrix(OpSum const &ops, Block const &block_in,
                             Block const &block_out, int64_t i0 = 0);
XDIAG_API csr_cx_mat csr_matrixC(Op const &op, Block const &block_in,
                                 Block const &block_out, int64_t i0 = 0);
XDIAG_API csr_cx_mat csr_matrixC(OpSum const &ops, Block const &block_in,
                                 Block const &block_out, int64_t i0 = 0);

template <typename idx_t, typename coeff_t>
XDIAG_API CSRMatrix<idx_t, coeff_t> csr_matrix(Op const &op, Block const &block,
                                               idx_t i0 = 0);
template <typename idx_t, typename coeff_t>
XDIAG_API CSRMatrix<idx_t, coeff_t>
csr_matrix(OpSum const &ops, Block const &block, idx_t i0 = 0);
template <typename idx_t, typename coeff_t>
XDIAG_API CSRMatrix<idx_t, coeff_t>
csr_matrix(Op const &op, Block const &block_in, Block const &block_out,
           idx_t i0 = 0);
template <typename idx_t, typename coeff_t>
XDIAG_API CSRMatrix<idx_t, coeff_t>
csr_matrix(OpSum const &ops, Block const &block_in, Block const &block_out,
           idx_t i0 = 0);

template <typename idx_t, typename coeff_t>
XDIAG_API arma::Mat<coeff_t> to_dense(CSRMatrix<idx_t, coeff_t> const &csr_mat);
} // namespace xdiag
