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

template <typename coeff_t = complex>
State<coeff_t> random_state(Block const &block, uint64_t seed = 42);
StateReal random_state_real(Block const &block, uint64_t seed = 42);
StateCplx random_state_cplx(Block const &block, uint64_t seed = 42);
  
} // namespace hydra
