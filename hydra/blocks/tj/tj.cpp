#include "tj.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/up_down_hole.h>

namespace hydra {

using namespace combinatorics;

template <class bit_t>
tJ<bit_t>::tJ(int n_sites, int nup, int ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      n_holes_(n_sites - charge_), size_holes_(binomial(n_sites, n_holes_)),
      size_spins_(binomial(charge_, nup)), size_(size_holes_ * size_spins_),
      lintable_holes_(n_sites, n_holes_), lintable_spins_(charge_, nup) {
  utils::check_nup_ndn_tj(n_sites, nup, ndn, "tJ");
}

template <class bit_t> bool tJ<bit_t>::operator==(tJ<bit_t> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_);
}

template <class bit_t> bool tJ<bit_t>::operator!=(tJ<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class tJ<uint16_t>;
template class tJ<uint32_t>;
template class tJ<uint64_t>;

} // namespace hydra
