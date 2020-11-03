#include <cassert>
#include <math.h>

#include "basis_tj.h"

#include <hydra/combinatorics/binomial.h>

namespace hydra {

template <class bit_t>
BasisTJ<bit_t>::BasisTJ(int const &n_sites, qn_t const &qn)
    : n_sites_(n_sites), qn_(qn) {
  assert(valid(qn, n_sites));
  BasisSpinHalf<bit_t> up(n_sites, qn.n_up);
  BasisSpinHalf<bit_t> holes(n_sites - qn.n_up, qn.n_dn);
  begin_ = BasisTJIterator<bit_t>(
      up, holes, {(*up.begin()).spins, (*holes.begin()).spins});
  end_ = BasisTJIterator<bit_t>(
      up, holes, {(*up.end()).spins, (*holes.begin()).spins});
}

template <class bit_t> int BasisTJ<bit_t>::n_sites() const {
  return n_sites_;
}

template <class bit_t> qn_tj BasisTJ<bit_t>::qn() const {
  return qn_;
}

template <class bit_t>
BasisTJIterator<bit_t> BasisTJ<bit_t>::begin() const {
  return begin_;
}

template <class bit_t>
BasisTJIterator<bit_t> BasisTJ<bit_t>::end() const {
  return end_;
}

template <class bit_t> uint64 BasisTJ<bit_t>::size() const {
  return combinatorics::binomial(n_sites_, qn_.n_up) *
    combinatorics::binomial(n_sites_ - qn_.n_up, qn_.n_dn);
}

template <class bit_t> uint64 BasisTJ<bit_t>::rawsize() const {
  return pow(3, n_sites_);
}

template <class bit_t>
BasisTJIterator<bit_t>::BasisTJIterator(
    BasisSpinHalf<bit_t> const &up, BasisSpinHalf<bit_t> const &holes,
    state_tj<bit_t> const &state)
    : n_sites_(up.n_sites()), holes_begin_(holes.begin()),
      holes_end_(holes.end()), holes_iter_(state_spinhalf<bit_t>({state.dns})),
      up_iter_(state_spinhalf<bit_t>({state.ups})) {
}

template class BasisTJ<uint16>;
template class BasisTJ<uint32>;
template class BasisTJ<uint64>;

template class BasisTJIterator<uint16>;
template class BasisTJIterator<uint32>;
template class BasisTJIterator<uint64>;

} // namespace hydra
