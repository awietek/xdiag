// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>
#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

#ifdef _OPENMP

template <typename idx_t, typename coeff_t, typename block_t>
void coo_matrix_fill(OpSum const &ops, block_t const &block_in,
                     block_t const &block_out,
                     std::vector<int64_t> const &nnz_thread, int64_t nnz,
                     idx_t *row, idx_t *col, coeff_t *data, idx_t i0);
#else

template <typename idx_t, typename coeff_t, typename block_t>
void coo_matrix_fill(OpSum const &ops, block_t const &block_in,
                     block_t const &block_out, int64_t nnz, idx_t *row,
                     idx_t *col, coeff_t *data, idx_t i0);
#endif

} // namespace xdiag::algebra
