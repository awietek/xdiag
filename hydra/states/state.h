#pragma once

#include <lila/all.h>

#include <hydra/common.h>

namespace hydra {

template <typename coeff_tt, class Block> class State {
public:
  using block_t = Block;
  using coeff_t = coeff_tt;

  State() = default;
  State(Block const &block, lila::Vector<coeff_t> const &vector)
      : block_(block), vector_(vector) {}
  State(lila::Vector<coeff_t> const &vector, Block const &block)
      : block_(block), vector_(vector) {}

  coeff_t operator()(idx_t idx) { return vector_(idx); }

  Block const &block() const { return block_; }
  lila::Vector<coeff_t> &vector() { return vector_; }
  lila::Vector<coeff_t> const &vector() const { return vector_; }
  idx_t size() const { return vector_.size(); }

private:
  Block block_;
  lila::Vector<coeff_t> vector_;
};

template <class Block> using StateReal = State<double, Block>;
template <class Block> using StateCplx = State<complex, Block>;

} // namespace hydra
