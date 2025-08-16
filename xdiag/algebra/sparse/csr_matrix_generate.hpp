// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>
#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
CSRMatrix<idx_t, coeff_t>
csr_matrix_generate(OpSum const &ops, block_t const &block_in,
                    block_t const &block_out, idx_t i0, bool transpose);

} // namespace xdiag::algebra
