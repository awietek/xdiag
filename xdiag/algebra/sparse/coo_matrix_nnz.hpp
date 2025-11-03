// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

#ifdef _OPENMP
template <typename coeff_t, typename block_t>
std::vector<int64_t> coo_matrix_nnz_thread(OpSum const &ops,
                                           block_t const &block_in,
                                           block_t const &block_out);

#else
template <typename coeff_t, typename block_t>
int64_t coo_matrix_nnz(OpSum const &ops, block_t const &block_in,
                       block_t const &block_out);

#endif

} // namespace xdiag::algebra
