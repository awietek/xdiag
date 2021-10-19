#include "tj.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/up_down_hole.h>
#include <hydra/models/utils/model_utils.h>

namespace hydra {

using namespace combinatorics;

template <class bit_t>
tJ<bit_t>::tJ(int n_sites, int nup, int ndn)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      size_up_(binomial(n_sites, nup)),
      size_holes_(binomial(n_sites - nup, ndn)), size_(size_up_ * size_holes_),
      dn_limits_for_up_(size_up_), dns_(size_), lintable_up_(n_sites, nup) {
  using combinatorics::Combinations;

  utils::check_nup_ndn_tj(n_sites, nup, ndn, "tJ");

  // Create dns array storing dn configs within the limits for upspin config
  idx_t idx = 0;
  idx_t idx_up = 0;
  for (bit_t up : Combinations<bit_t>(n_sites, nup)) {
    dn_limits_for_up_[idx_up].first = idx;
    for (bit_t holes : Combinations<bit_t>(n_sites - nup, ndn)) {
      bit_t dn = up_hole_to_down(up, holes);
      dns_[idx] = dn;
      ++idx;
    }
    dn_limits_for_up_[idx_up].second = idx;
    ++idx_up;
  }
  assert(idx_up == size_up_);
  assert(idx == size_);
}

template <class bit_t> idx_t tJ<bit_t>::index(bit_t up, bit_t dn) const {
  idx_t up_idx = lintable_up_.index(up);
  auto [dn_lower, dn_upper] = dn_limits_for_up_[up_idx];
  auto it =
      std::lower_bound(dns_.begin() + dn_lower, dns_.begin() + dn_upper, dn);
  if ((it == dns_.begin() + dn_upper) || (*it != dn))
    return invalid_index;
  else
    return std::distance(dns_.begin(), it);
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

template class tJ<uint16>;
template class tJ<uint32>;
template class tJ<uint64>;

} // namespace hydra
