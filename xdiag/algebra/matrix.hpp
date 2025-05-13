// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

XDIAG_API arma::mat matrix(Op const &op, Block const &block);
XDIAG_API arma::mat matrix(OpSum const &ops, Block const &block);
XDIAG_API arma::cx_mat matrixC(Op const &op, Block const &block);
XDIAG_API arma::cx_mat matrixC(OpSum const &ops, Block const &block);

XDIAG_API arma::mat matrix(Op const &op, Block const &block_in,
                           Block const &block_out);
XDIAG_API arma::mat matrix(OpSum const &ops, Block const &block_in,
                           Block const &block_out);
XDIAG_API arma::cx_mat matrixC(Op const &op, Block const &block_in,
                               Block const &block_out);
XDIAG_API arma::cx_mat matrixC(OpSum const &ops, Block const &block_in,
                               Block const &block_out);

// developer method
template <typename coeff_t, class block_t>
XDIAG_API void matrix(coeff_t *mat, OpSum const &ops, block_t const &block_in,
                      block_t const &block_out);

} // namespace xdiag
