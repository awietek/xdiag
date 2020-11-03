#include <cassert>
#include <math.h>
#include <sstream>

#include "basis_spinhalf.h"

#include <hydra/utils/bitops.h>
#include <hydra/combinatorics/binomial.h>

namespace hydra {

template <class bit_t>
BasisSpinHalf<bit_t>::BasisSpinHalf(int const &n_sites, qn_t const &qn)
    : n_sites_(n_sites), qn_(qn) {
  assert(valid(qn, n_sites_));
  bit_t begin_state = (((bit_t)1 << qn.n_up) - 1);
  bit_t end_state = begin_state << (n_sites_ - qn.n_up);
  end_state = combinatorics::get_next_pattern<bit_t>(end_state);
  begin_ = BasisSpinHalfIterator<bit_t>(state_t({begin_state}));
  end_ = BasisSpinHalfIterator<bit_t>(state_t({end_state}));
}

template <class bit_t>
BasisSpinHalf<bit_t>::BasisSpinHalf(int const &n_sites, int const &n_up)
    : BasisSpinHalf(n_sites, qn_t({n_up})) {}

template <class bit_t> int BasisSpinHalf<bit_t>::n_sites() const {
  return n_sites_;
}

template <class bit_t> qn_spinhalf BasisSpinHalf<bit_t>::qn() const {
  return qn_;
}

template <class bit_t>
BasisSpinHalfIterator<bit_t> BasisSpinHalf<bit_t>::begin() const {
  return begin_;
}

template <class bit_t>
BasisSpinHalfIterator<bit_t> BasisSpinHalf<bit_t>::end() const {
  return end_;
}

template <class bit_t> uint64 BasisSpinHalf<bit_t>::size() const {
  return combinatorics::binomial(n_sites_, qn_.n_up);
}

template <class bit_t> uint64 BasisSpinHalf<bit_t>::rawsize() const {
  return pow(2, n_sites_);
}

template <class bit_t>
BasisSpinHalfIterator<bit_t>::BasisSpinHalfIterator(const state_t &state)
    : current_(state.spins) {}

template class BasisSpinHalf<uint16>;
template class BasisSpinHalf<uint32>;
template class BasisSpinHalf<uint64>;

template class BasisSpinHalfIterator<uint16>;
template class BasisSpinHalfIterator<uint32>;
template class BasisSpinHalfIterator<uint64>;

} // namespace hydra
