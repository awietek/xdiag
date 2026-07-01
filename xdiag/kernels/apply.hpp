// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Computes vec_out = ops * vec_in where vec_in/vec_out are arma::Col or
// arma::Mat (real or complex).
// vec_out is zeroed before accumulation. block_in and block_out must be the
// same Block type.
template <typename mat_t>
XDIAG_API void apply(Op const &op, Block const &block_in, mat_t const &vec_in,
                     Block const &block_out, mat_t &vec_out);
template <typename mat_t>
XDIAG_API void apply(Monomial const &mono, Block const &block_in,
                     mat_t const &vec_in, Block const &block_out,
                     mat_t &vec_out);
template <typename mat_t>
XDIAG_API void apply(OpSum const &ops, Block const &block_in,
                     mat_t const &vec_in, Block const &block_out,
                     mat_t &vec_out);

} // namespace xdiag
