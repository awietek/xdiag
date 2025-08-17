// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

#include <xdiag/blocks/blocks.hpp>

namespace xdiag {

XDIAG_API State apply(Op const &op, State const &v);
XDIAG_API State apply(OpSum const &ops, State const &v);
XDIAG_API void apply(Op const &op, State const &v, State &w);
XDIAG_API void apply(OpSum const &ops, State const &v, State &w);

XDIAG_API void apply(OpSum const &ops, Block const &block_in,
                     arma::vec const &vec_in, Block const &block_out,
                     arma::vec &vec_out);
XDIAG_API void apply(OpSum const &ops, Block const &block_in,
                     arma::cx_vec const &vec_in, Block const &block_out,
                     arma::cx_vec &vec_out);
XDIAG_API void apply(OpSum const &ops, Block const &block_in,
                     arma::mat const &mat_in, Block const &block_out,
                     arma::mat &mat_out);
XDIAG_API void apply(OpSum const &ops, Block const &block_in,
                     arma::cx_mat const &mat_in, Block const &block_out,
                     arma::cx_mat &mat_out);
  
template <typename mat_t, typename block_t>
XDIAG_API void apply(OpSum const &ops, block_t const &block_in,
                     mat_t const &mat_in, block_t const &block_out,
                     mat_t &mat_out);

} // namespace xdiag
