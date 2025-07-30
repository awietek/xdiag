// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
COOMatrix<idx_t, coeff_t>
coo_matrix_generate(OpSum const &ops, block_t const &block_in,
                    block_t const &block_out, idx_t i0);

} // namespace xdiag::algebra
