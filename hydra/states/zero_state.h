#pragma once

#include "extern/armadillo/armadillo"

namespace hydra {

template <class State> void ZeroState(State &state) {
  state.vector().zeros();
}

template <class coeff_t = complex, class Block>
State<coeff_t, Block> ZeroState(Block const &block) {
  auto v = arma::zeros<arma::Col<coeff_t>>(block.size());
  return State(block, v);
}

template <class Block> StateReal<Block> ZeroStateReal(Block const &block) {
  return ZeroState<double, Block>(block);
}

template <class Block> StateCplx<Block> ZeroStateCplx(Block const &block) {
  return ZeroState<complex, Block>(block);
}

} // namespace hydra
