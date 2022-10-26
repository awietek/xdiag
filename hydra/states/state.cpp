#include "state.h"

namespace hydra {

template <typename coeff_tt>
State<coeff_t>::State(Block const &block)
    : block_(block), vector_(size(block), fill::zeros) {}

template <typename coeff_tt>
State<coeff_t>::State(Block const &block, arma::Col<coeff_t> const &vector)
    : block_(block), vector_(vector) {}

template <typename coeff_tt>
State<coeff_t>::State(Block const &block, ProductState const &pstate)
    : block_(block), vector_(size(block), fill::zeros) {
  fill(pstate, *this);
}

template <typename coeff_tt>
State<coeff_t>::State(Block const &block, RandomState const &rstate)
    : block_(block), vector_(size(block), fill::zeros) {
  fill(pstate, *this);
}

template class State<double>;
template class State<complex>;
  
} // namespace hydra
