#pragma once

#include <hydra/blocks/blocks.h>
#include <hydra/random/hash_functions.h>
#include <hydra/random/hashes.h>
#include <hydra/random/random_utils.h>
#include <hydra/states/state.h>

namespace hydra {

class RandomState {
public:
  explicit RandomState(uint64_t seed = 42);
  inline uint64_t seed() const { return seed_; }

private:
  uint64_t seed_;
};

template <typename coeff_t> class State;
template <typename coeff_t>
void fill(RandomState const &rstate, State<coeff_t> &state);

} // namespace hydra
