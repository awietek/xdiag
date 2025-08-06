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

// int64_t, double (default)
template <typename idx_t = int64_t, typename coeff_t = double>
XDIAG_API CSCMatrix<idx_t, coeff_t>
csc_matrix(OpSum const &ops, Block const &block, idx_t i0 = 0);

template <typename idx_t = int64_t, typename coeff_t = double>
XDIAG_API CSCMatrix<idx_t, coeff_t>
csc_matrix(OpSum const &ops, Block const &block_in, Block const &block_out,
           idx_t i0 = 0);

// int64_t, complex
XDIAG_API CSCMatrix<int64_t, complex>
csc_matrixC(OpSum const &ops, Block const &block, int64_t i0 = 0);
XDIAG_API CSCMatrix<int64_t, complex> csc_matrixC(OpSum const &ops,
                                                  Block const &block_in,
                                                  Block const &block_out,
                                                  int64_t i0 = 0);
// int32_t, double
XDIAG_API CSCMatrix<int32_t, double>
csc_matrix_32(OpSum const &ops, Block const &block, int32_t i0 = 0);
XDIAG_API CSCMatrix<int32_t, double> csc_matrix_32(OpSum const &ops,
                                                   Block const &block_in,
                                                   Block const &block_out,
                                                   int32_t i0 = 0);
// int32_t, complex
XDIAG_API CSCMatrix<int32_t, complex>
csc_matrixC_32(OpSum const &ops, Block const &block, int32_t i0 = 0);
XDIAG_API CSCMatrix<int32_t, complex> csc_matrixC_32(OpSum const &ops,
                                                     Block const &block_in,
                                                     Block const &block_out,
                                                     int32_t i0 = 0);

template <typename idx_t, typename coeff_t>
XDIAG_API arma::Mat<coeff_t> to_dense(CSCMatrix<idx_t, coeff_t> const &csc_mat);
} // namespace xdiag
