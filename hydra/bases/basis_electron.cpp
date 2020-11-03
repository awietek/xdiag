#include <cassert>
#include <math.h>
#include <string>

#include "basis_electron.h"

#include <hydra/combinatorics/binomial.h>

namespace hydra {

template <class bit_t>
BasisElectron<bit_t>::BasisElectron(int const &n_sites, qn_t const &qn)
    : n_sites_(n_sites), qn_(qn) {
  assert(valid(qn, n_sites));
  BasisSpinHalf<bit_t> up(n_sites, qn.n_up);
  BasisSpinHalf<bit_t> down(n_sites, qn.n_dn);
  begin_ = BasisElectronIterator<bit_t>(
      up, down, {(*up.begin()).spins, (*down.begin()).spins});
  end_ = BasisElectronIterator<bit_t>(
      up, down, {(*up.end()).spins, (*down.begin()).spins});
}

template <class bit_t> int BasisElectron<bit_t>::n_sites() const {
  return n_sites_;
}

template <class bit_t> qn_electron BasisElectron<bit_t>::qn() const {
  return qn_;
}

template <class bit_t>
BasisElectronIterator<bit_t> BasisElectron<bit_t>::begin() const {
  return begin_;
}

template <class bit_t>
BasisElectronIterator<bit_t> BasisElectron<bit_t>::end() const {
  return end_;
}

template <class bit_t> uint64 BasisElectron<bit_t>::size() const {
  return combinatorics::binomial(n_sites_, qn_.n_up) *
         combinatorics::binomial(n_sites_, qn_.n_dn);
}

template <class bit_t> uint64 BasisElectron<bit_t>::rawsize() const {
  return pow(4, n_sites_);
}

template <class bit_t>
BasisElectronIterator<bit_t>::BasisElectronIterator(
    BasisSpinHalf<bit_t> const &up, BasisSpinHalf<bit_t> const &down,
    state_electron<bit_t> const &state)
    : n_sites_(up.n_sites()), down_begin_(down.begin()), down_end_(down.end()),
      down_iter_(state_spinhalf<bit_t>({state.dns})),
      up_iter_(state_spinhalf<bit_t>({state.ups})) {
  assert(up.n_sites() == down.n_sites());
}

template class BasisElectron<uint16>;
template class BasisElectron<uint32>;
template class BasisElectron<uint64>;

template class BasisElectronIterator<uint16>;
template class BasisElectronIterator<uint32>;
template class BasisElectronIterator<uint64>;

} // namespace hydra
