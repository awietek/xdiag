#pragma once

#include <hydra/blocks/blocks.h>
#include <hydra/random/hash_functions.h>
#include <hydra/random/hashes.h>
#include <hydra/random/random_utils.h>
#include <hydra/states/state.h>

namespace hydra {

template <class coeff_t, class Block>
inline void random_state(State<coeff_t, Block> &state, uint32_t seed = 42) {
  if constexpr (mpi::is_mpi_block<Block>()) {
    seed += 0x01000193 * state.block.mpi_rank();
  }
  uint32_t seed_modified =
      random::hash_combine(seed, random::hash(state.block()));
  random::fill_random_normal_vector(state.vector(), seed_modified);
  coeff_t norm = arma::norm(state.vector());
  state.vector() /= norm;
}

template <class coeff_t = complex, class Block>
inline State<coeff_t, Block> random_state(Block const &block,
                                          uint32_t seed = 42) {
  if constexpr (mpi::is_mpi_block<Block>) {
    seed += 0x01000193 * block.mpi_rank();
  }

  arma::Col<coeff_t> v(block.size(), arma::fill::zeros);
  uint32_t seed_modified = random::hash_combine(seed, random::hash(block));
  random::fill_random_normal_vector(v, seed_modified);
  coeff_t norm = arma::norm(v);
  v /= norm;
  return State(block, v);
}
template <class Block>
inline StateReal<Block> random_state_real(Block const &block,
                                          uint32_t seed = 42) {
  return random_state<double, Block>(block, seed);
}

template <class Block>
inline StateCplx<Block> random_state_cplx(Block const &block,
                                          uint32_t seed = 42) {
  return random_state<complex, Block>(block, seed);
}

} // namespace hydra
