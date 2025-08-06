// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>
#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
std::vector<idx_t> csr_matrix_nnz(OpSum const &ops, block_t const &block_in,
                                  block_t const &block_out, bool transpose);

} // namespace xdiag::algebra
