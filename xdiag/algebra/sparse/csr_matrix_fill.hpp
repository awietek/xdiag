// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>
#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
void csr_matrix_fill(OpSum const &ops, block_t const &block_in,
                     block_t const &block_out, int64_t nnz,
                     std::vector<idx_t> const &n_elements_in_row, idx_t *rowptr,
                     idx_t *col, coeff_t *data, idx_t i0, bool transpose);

} // namespace xdiag::algebra
