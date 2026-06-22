// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Computes vec_out = ops * vec_in, where ops is an Op, Monomial, or OpSum,
// and vec_in/vec_out are arma::Col or arma::Mat (real or complex).
// vec_out is zeroed before accumulation. block_in and block_out must be the
// same Block type; for hermitian operators acting within one sector they are
// identical.
template <typename op_t, typename mat_t>
XDIAG_API void apply(op_t const &ops, Block const &block_in,
                     mat_t const &vec_in, Block const &block_out,
                     mat_t &vec_out);

} // namespace xdiag
