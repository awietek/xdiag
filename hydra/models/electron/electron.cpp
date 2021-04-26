#include "electron.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/models/electron_utils.h>

namespace hydra {

template <class bit_t>
Electron<bit_t>::Electron(int n_sites, int nup, int ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      lintable_up_(n_sites, n_up_), lintable_dn_(n_sites, n_dn_),
      size_up_(combinatorics::binomial(n_sites, n_up_)),
      size_dn_(combinatorics::binomial(n_sites, n_dn_)),
      size_(size_up_ * size_dn_) {
  electron::check_nup_ndn(n_sites, nup, ndn);
}

template <class bit_t>
bool Electron<bit_t>::operator==(Electron<bit_t> const &rhs) const{
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_);
}

template <class bit_t>
bool Electron<bit_t>::operator!=(Electron<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class Electron<uint16>;
template class Electron<uint32>;
template class Electron<uint64>;

} // namespace hydra
