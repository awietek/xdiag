#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

State product(Block const &block, std::vector<std::string> const &local_state,
              bool real = true);
template <typename block_t>
State product(block_t const &block, std::vector<std::string> const &local_state,
              bool real = true);

State rand(Block const &block, bool real = true, int64_t seed = 42,
           bool normalized = true);
template <typename block_t>
State rand(block_t const &block, bool real = true, int64_t seed = 42,
           bool normalized = true);

State zeros(Block const &block, bool real = true, int64_t n_cols = 1);
template <typename block_t>
State zeros(block_t const &block, bool real = true, int64_t n_cols = 1);
void zero(State &state);

} // namespace xdiag
