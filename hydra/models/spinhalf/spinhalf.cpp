#include "spinhalf.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/models/utils/model_utils.h>

namespace hydra {

template <class bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites, int n_up)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_), lintable_(n_sites, n_up),
      size_(combinatorics::binomial(n_sites, n_up)) {
  utils::check_nup_spinhalf(n_sites, n_up, "Spinhalf");
}

template <class bit_t>
bool Spinhalf<bit_t>::operator==(Spinhalf<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_);
}

template <class bit_t>
bool Spinhalf<bit_t>::operator!=(Spinhalf<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class Spinhalf<uint16>;
template class Spinhalf<uint32>;
template class Spinhalf<uint64>;

} // namespace hydra
