#pragma once

#include <lila/all.h>

namespace hydra {

template <typename coeff_t, class Block> class State {
public:
  State() = default;
  State(Block const &block, lila::Vector<coeff_t> const &vector)
      : block_(block), vector_(vector) {}
  State(lila::Vector<coeff_t> const &vector, Block const &block)
      : block_(block), vector_(vector) {}

  Block const &block() const { return block_; }
  lila::Vector<coeff_t> &vector() { return vector_; }
  lila::Vector<coeff_t> const &vector() const { return vector_; }

private:
  Block block_;
  lila::Vector<coeff_t> vector_;
};

template <class Block> using StateReal = State<double, Block>;
template <class Block> using StateCplx = State<complex, Block>;

template <class coeff_t, class Block>
State<coeff_t, Block> ZeroState(Block const &block) {
  auto v = lila::Zeros<coeff_t>(block.size());
  return State(block, v);
}

template <class Block> StateReal<Block> ZeroStateReal(Block const &block) {
  return ZeroState<double, Block>(block);
}

template <class Block> StateCplx<Block> ZeroStateCplx(Block const &block) {
  return ZeroState<complex, Block>(block);
}

template <class coeff_t, class Block>
State<coeff_t, Block> RandomState(Block const &block, int seed = 42) {
  if constexpr (detail::is_mpi_block<Block>) {
    seed += 0x01000193 * block.mpi_rank();
  }
  auto v = lila::Random<coeff_t>(block.size());
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
