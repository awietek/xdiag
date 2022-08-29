#pragma once

#include <lila/all.h>

#include <hydra/blocks/blocks.h>
#include <hydra/states/state.h>
#include <hydra/utils/random_utils.h>

namespace hydra {

template <class State> void RandomState(State &state, int seed = 42) {
  if constexpr (detail::is_mpi_block<typename State::block_t>()) {
    seed += 0x01000193 * state.block.mpi_rank();
  }
  random::fill_random_normal_vector(state.vector(), seed);
}

template <class coeff_t, class Block>
State<coeff_t, Block> RandomState(Block const &block, int seed = 42) {
  if constexpr (detail::is_mpi_block<Block>) {
    seed += 0x01000193 * block.mpi_rank();
  }

  auto v = lila::Zeros<coeff_t>(block.size());
  random::fill_random_normal_vector(v, seed);
  return State(block, v);
}
template <class Block>
StateReal<Block> RandomStateReal(Block const &block, int seed = 42) {
  return RandomState<double, Block>(block);
}

template <class Block>
StateCplx<Block> RandomStateCplx(Block const &block, int seed = 42) {
  return RandomState<complex, Block>(block);
}

} // namespace hydra
