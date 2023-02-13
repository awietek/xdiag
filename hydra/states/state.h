#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>

namespace hydra {

class ProductState;
class RandomState;

template <typename coeff_tt = complex> class State {
public:
  using coeff_t = coeff_tt;

  State() = default;
  State(Block const &block);
  State(Block const &block, arma::Col<coeff_t> const &vector);
  State(Block const &block, ProductState const &pstate);
  State(Block const &block, RandomState const &rstate);

  coeff_t operator()(idx_t idx) { return vector_(idx); }

  Block const &block() const { return block_; }
  arma::Col<coeff_t> &vector() { return vector_; }
  arma::Col<coeff_t> const &vector() const { return vector_; }
  idx_t size() const { return vector_.size(); }

private:
  Block block_;
  arma::Col<coeff_t> vector_;
};

using StateReal = State<double>;
using StateCplx = State<complex>;

StateCplx to_cplx(StateReal const &state);

} // namespace hydra
