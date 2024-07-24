#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

class RandomState {
public:
  RandomState(int64_t seed = 42, bool normalized = true);
  int64_t seed() const;
  bool normalized() const;

private:
  int64_t seed_;
  bool normalized_;
};

void fill(State &state, RandomState const &rstate, int64_t col = 0);

State rand(Block const &block, bool real = true, int64_t seed = 42,
           bool normalized = true);
template <typename block_t>
State rand(block_t const &block, bool real = true, int64_t seed = 42,
           bool normalized = true);

std::ostream &operator<<(std::ostream &out, RandomState const &state);
std::string to_string(RandomState const &state);

} // namespace xdiag
