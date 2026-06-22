// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>
#include <xdiag/blocks/blocks.hpp>

namespace xdiag {

// Validates the inputs of a sparse-matrix construction before building:
//   - the zero index i0 is 0 or 1,
//   - ops maps block_in into block_out's symmetry sector (blocks_match),
//   - a real matrix (coeff_t real) is only requested for real ops and blocks,
//   - the block dimensions fit into the index type idx_t.
// Throws on any violation. Defined for each (idx_t, coeff_t, block_t) used by
// the sparse-matrix builders (see explicit instantiations in valid.cpp).
template <typename idx_t, typename coeff_t>
void check_valid_sparse_matrix(OpSum const &ops, Block const &block_in,
                               Block const &block_out, idx_t i0);

} // namespace xdiag
