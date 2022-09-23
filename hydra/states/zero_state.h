#pragma once

#include "extern/armadillo/armadillo"

namespace hydra {

template <class State> void zero_state(State &state) { state.vector().zeros(); }

template <class coeff_t = complex, class Block>
State<coeff_t, Block> zero_state(Block const &block) {
  arma::Col<coeff_t> v(block.size(), arma::fill::zeros);
  return State(block, v);
}

template <class Block> StateReal<Block> zero_state_real(Block const &block) {
  return zero_state<double, Block>(block);
}

template <class Block> StateCplx<Block> ZeroStateCplx(Block const &block) {
  return zero_state<complex, Block>(block);
}

} // namespace hydra
