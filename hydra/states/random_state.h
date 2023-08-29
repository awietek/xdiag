#pragma once

#include <hydra/common.h>
#include <hydra/states/state.h>
#include <hydra/blocks/blocks.h>

namespace hydra {

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

State random_state(block_variant_t const &block, bool real = true,
                   int64_t seed = 42, bool normalized = true);
template <typename block_t>
State random_state(block_t const &block, bool real = true, int64_t seed = 42,
                   bool normalized = true);

} // namespace hydra
