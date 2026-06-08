// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

// Build a State holding the single product state given by its per-site local
// quantum-number indices (spin-1/2: Dn=0, Up=1; boson: occupation number).
XDIAG_API State product_state(Block const &block,
                              std::vector<int64_t> const &local_state,
                              bool real = true);
template <typename block_t>
XDIAG_API State product_state(block_t const &block,
                              std::vector<int64_t> const &local_state,
                              bool real = true);

XDIAG_API State random_state(Block const &block, bool real = true,
                             int64_t ncols = 1, int64_t seed = 42,
                             bool normalized = true);
template <typename block_t>
XDIAG_API State random_state(block_t const &block, bool real = true,
                             int64_t ncols = 1, int64_t seed = 42,
                             bool normalized = true);

XDIAG_API State zero_state(Block const &block, bool real = true,
                           int64_t ncols = 1);
template <typename block_t>
XDIAG_API State zero_state(block_t const &block, bool real = true,
                           int64_t ncols = 1);
XDIAG_API void zero(State &state);

} // namespace xdiag
