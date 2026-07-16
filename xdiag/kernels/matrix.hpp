// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Returns the dense matrix representation of ops acting from block_in to
// block_out. The output block is inferred from the operator quantum numbers.
// op_t may be Op, Monomial, or OpSum. Use matrixC to force a complex result.
template <typename op_t>
XDIAG_API arma::mat matrix(op_t const &op, Block const &block);

template <typename op_t>
XDIAG_API arma::cx_mat matrixC(op_t const &op, Block const &block);

// Two-block overloads: explicitly specify both input and output block.
// Useful for off-diagonal sectors (e.g. hopping between particle-number
// sectors) where the output block differs from the input block.
template <typename op_t>
XDIAG_API arma::mat matrix(op_t const &op, Block const &block_in,
                           Block const &block_out);
template <typename op_t>
XDIAG_API arma::cx_mat matrixC(op_t const &op, Block const &block_in,
                               Block const &block_out);

// Developer overload used by the Julia wrapper: fills a caller-allocated,
// column-major buffer mat of size size(block_out) x dim(block_in). The buffer
// is zeroed by this routine before the kernel accumulates into it (the kernel
// uses +=), so callers may pass uninitialized memory. mat must stay valid for
// the duration of the call only; no pointer to it is retained. coeff_t must be
// double or complex. Throws for distributed blocks (no dense representation).
template <typename coeff_t>
XDIAG_API void matrix(OpSum const &ops, Block const &block_in,
                      Block const &block_out, coeff_t *mat);

} // namespace xdiag
