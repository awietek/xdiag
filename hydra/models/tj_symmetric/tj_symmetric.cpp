#include "electron_symmetric.h"

#include <cassert>
#include <tuple>

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>

#include <hydra/models/electron/electron_utils.h>

namespace hydra {
namespace electron {

template <class bit_t, class SymmetryGroup>
double compute_norm(bit_t ups, bit_t dns, SymmetryGroup &&symmetry_group,
                    Representation const &irrep) {
  assert(symmetry_group.size() == irrep.size());
  complex amplitude = 0.0;
  for (int sym = 0; sym < (int)symmetry_group.size(); ++sym) {
    bit_t tups = symmetry_group.apply(sym, ups);
    if (tups == ups) {
      bit_t tdns = symmetry_group.apply(sym, dns);
      double fermi_sign_ups = symmetry_group.fermi_sign(sym, ups);
      if (tdns == dns) {
        double fermi_sign_dns = symmetry_group.fermi_sign(sym, dns);
        amplitude += fermi_sign_ups * fermi_sign_dns * irrep.character(sym);
      }
    }
  }
  return std::sqrt(std::abs(amplitude));
}

} // namespace electron

template <class bit_t, class SymmetryGroup>
tJSymmetric<bit_t, SymmetryGroup>::tJSymmetric(
    int n_sites, int nup, int ndn, SymmetryGroup symmetry_group,
    Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      symmetry_group_(symmetry_group), irrep_(irrep) {

  electron::check_nup_ndn(n_sites, nup, ndn);

  if (irrep.allowed_symmetries().size() > 0) {
    symmetry_group_ = symmetry_group.subgroup(irrep.allowed_symmetries());
  }
  
  // Compute downspin configurations and up limits
  idx_t idx = 0;
  for (auto ups : Combinations(n_sites, nup)) {

    if (symmetry_group_.representative(ups) != ups) continue;
    
    idx_t lower = idx;
    for (auto dns : Combinations(n_sites, ndn)) {

      // Determine representative in up/dn ordering
      auto [rep_ups, rep_dns] = representative(ups, dns);

      // If state is a representative ...
      if ((rep_ups == ups) && (rep_dns == dns)) {
        double norm =
            electron::compute_norm<bit_t>(ups, dns, symmetry_group_, irrep);

        // ... and norm is nonzero, register the state and its norm
        if (norm > 1e-6) {  // tolerance big as 1e-6 since root is taken
          dns_.push_back(rep_dns);
          norms_.push_back(norm);
          ++idx;
        }
      }
    }
    idx_t upper = idx;
    if (upper > lower)
      ups_lower_upper_[ups] = {lower, upper};
  }

  // Compute upspin configurations and dn limits
  idx = 0;
  for (auto dns : Combinations(n_sites, ndn)) {

    if (symmetry_group_.representative(dns) != dns) continue;
    
    idx_t lower = idx;
    for (auto ups : Combinations(n_sites, nup)) {

      // Determine representative in dn/up ordering (switch)
      auto [rep_dns_switch, rep_ups_switch] = representative(dns, ups);

      // If state is a (switch) representative ...
      if ((rep_ups_switch == ups) && (rep_dns_switch == dns)) {
        double norm =
            electron::compute_norm<bit_t>(ups, dns, symmetry_group_, irrep);

        // ... and has non-zero norm
        if (std::abs(norm) > 1e-6) {   // tolerance big as 1e-6 since root is taken

          // keep the up configuration
          ups_.push_back(rep_ups_switch);

          // compute its actual representative
          auto [rep_ups, rep_dns, sym] =
              representative_index(rep_ups_switch, rep_dns_switch);
          index_switch_to_index_.push_back(index(rep_ups, rep_dns));
          complex chi = irrep.character(sym) *
              symmetry_group_.fermi_sign(sym, rep_ups_switch) *
              symmetry_group_.fermi_sign(sym, rep_dns_switch);
          character_switch_.push_back(chi);
          ++idx;
        }
      }
    }
    idx_t upper = idx;
    if (upper > lower)
      dns_lower_upper_[dns] = {lower, upper};
  }
  assert(ups_.size() == dns_.size());
  size_ = (idx_t)ups_.size();
}

template <class bit_t, class SymmetryGroup>
std::tuple<bit_t, bit_t>
tJSymmetric<bit_t, SymmetryGroup>::representative(bit_t ups,
                                                        bit_t dns) const {
  auto [rep_ups, n_sym, sym_ptr] = symmetry_group_.representative_indices(ups);
  bit_t rep_dns = std::numeric_limits<bit_t>::max();
  for (int i = 0; i < n_sym; ++i) {
    bit_t tdns = symmetry_group_.apply(sym_ptr[i], dns);
    if (tdns < rep_dns)
      rep_dns = tdns;
  }
  return {rep_ups, rep_dns};
}

template <class bit_t, class SymmetryGroup>
std::tuple<bit_t, bit_t, int>
tJSymmetric<bit_t, SymmetryGroup>::representative_index(bit_t ups,
                                                              bit_t dns) const {
  auto [rep_ups, n_sym, sym_ptr] = symmetry_group_.representative_indices(ups);
  bit_t rep_dns = std::numeric_limits<bit_t>::max();
  int sym = 0;
  for (int i = 0; i < n_sym; ++i) {
    bit_t tdns = symmetry_group_.apply(sym_ptr[i], dns);
    if (tdns < rep_dns) {
      rep_dns = tdns;
      sym = sym_ptr[i];
    }
  }
  return {rep_ups, rep_dns, sym};
}

template <class bit_t, class SymmetryGroup>
idx_t tJSymmetric<bit_t, SymmetryGroup>::index(bit_t ups,
                                                     bit_t dns) const {
  auto it1 = ups_lower_upper_.find(ups);
  if (it1 == ups_lower_upper_.end())
    return invalid_index;
  auto [lower, upper] = it1->second;
  auto it = std::lower_bound(dns_.begin() + lower, dns_.begin() + upper, dns);
  if ((it == dns_.begin() + upper) || (*it != dns))
    return invalid_index;
  else
    return std::distance(dns_.begin(), it);
}

template <class bit_t, class SymmetryGroup>
idx_t tJSymmetric<bit_t, SymmetryGroup>::index_switch(bit_t ups,
                                                            bit_t dns) const {
  auto it1 = dns_lower_upper_.find(dns);
  if (it1 == dns_lower_upper_.end())
    return invalid_index;
  auto [lower, upper] = it1->second;
  auto it = std::lower_bound(ups_.begin() + lower, ups_.begin() + upper, ups);
  if ((it == ups_.begin() + upper) || (*it != ups))
    return invalid_index;
  else
    return std::distance(ups_.begin(), it);
}

template <class bit_t, class SymmetryGroup>
bool tJSymmetric<bit_t, SymmetryGroup>::operator==(
    tJSymmetric<bit_t, SymmetryGroup> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (symmetry_group_ == rhs.symmetry_group_) && (irrep_ == rhs.irrep_);
}

template <class bit_t, class SymmetryGroup>
bool tJSymmetric<bit_t, SymmetryGroup>::operator!=(
    tJSymmetric<bit_t, SymmetryGroup> const &rhs) const {
  return !operator==(rhs);
}

template class tJSymmetric<uint16, SpaceGroup<uint16>>;
template class tJSymmetric<uint32, SpaceGroup<uint32>>;
template class tJSymmetric<uint64, SpaceGroup<uint64>>;

} // namespace hydra
