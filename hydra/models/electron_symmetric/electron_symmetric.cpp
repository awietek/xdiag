#include "electron_symmetric.h"

#include <cassert>
#include <tuple>

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>

#include <hydra/models/utils/model_utils.h>
#include <hydra/models/utils/symmetrized_norm.h>

namespace hydra {

template <class bit_t, class GroupAction>
ElectronSymmetric<bit_t, GroupAction>::ElectronSymmetric(
    int n_sites, int nup, int ndn, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      permutation_group_(permutation_group), irrep_(irrep) {

  utils::check_nup_ndn_electron(n_sites, nup, ndn, "ElectronSymmetric");

  if (irrep.allowed_symmetries().size() > 0) {
    permutation_group_ = permutation_group.subgroup(irrep.allowed_symmetries());
  }
  group_action_ = GroupAction(permutation_group_);

  // Compute downspin configurations and up limits
  idx_t idx = 0;
  for (auto ups : Combinations(n_sites, nup)) {

    if (group_action_.representative(ups) != ups)
      continue;

    idx_t lower = idx;
    for (auto dns : Combinations(n_sites, ndn)) {

      // Determine representative in up/dn ordering
      auto [rep_ups, rep_dns] = representative(ups, dns);

      // If state is a representative ...
      if ((rep_ups == ups) && (rep_dns == dns)) {
        double norm =
            utils::symmetrized_norm_electron(ups, dns, group_action_, irrep);

        // ... and norm is nonzero, register the state and its norm
        if (norm > 1e-6) { // tolerance big as 1e-6 since root is taken
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

    if (group_action_.representative(dns) != dns)
      continue;

    idx_t lower = idx;
    for (auto ups : Combinations(n_sites, nup)) {

      // Determine representative in dn/up ordering (switch)
      auto [rep_dns_switch, rep_ups_switch] = representative(dns, ups);

      // If state is a (switch) representative ...
      if ((rep_ups_switch == ups) && (rep_dns_switch == dns)) {
        double norm =
            utils::symmetrized_norm_electron(ups, dns, group_action_, irrep);

        // ... and has non-zero norm
        if (std::abs(norm) >
            1e-6) { // tolerance big as 1e-6 since root is taken

          // keep the up configuration
          ups_.push_back(rep_ups_switch);

          // compute its actual representative
          auto [rep_ups, rep_dns, sym] =
              representative_index(rep_ups_switch, rep_dns_switch);
          index_switch_to_index_.push_back(index(rep_ups, rep_dns));
          complex chi = irrep.character(sym) *
                        group_action_.fermi_sign(sym, rep_ups_switch) *
                        group_action_.fermi_sign(sym, rep_dns_switch);
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

template <class bit_t, class GroupAction>
std::tuple<bit_t, bit_t>
ElectronSymmetric<bit_t, GroupAction>::representative(bit_t ups,
                                                      bit_t dns) const {
  auto [rep_ups, n_sym, sym_ptr] = group_action_.representative_indices(ups);
  bit_t rep_dns = std::numeric_limits<bit_t>::max();
  for (int i = 0; i < n_sym; ++i) {
    bit_t tdns = group_action_.apply(sym_ptr[i], dns);
    if (tdns < rep_dns)
      rep_dns = tdns;
  }
  return {rep_ups, rep_dns};
}

template <class bit_t, class GroupAction>
std::tuple<bit_t, bit_t, int>
ElectronSymmetric<bit_t, GroupAction>::representative_index(bit_t ups,
                                                            bit_t dns) const {
  auto [rep_ups, n_sym, sym_ptr] = group_action_.representative_indices(ups);
  bit_t rep_dns = std::numeric_limits<bit_t>::max();
  int sym = 0;
  for (int i = 0; i < n_sym; ++i) {
    bit_t tdns = group_action_.apply(sym_ptr[i], dns);
    if (tdns < rep_dns) {
      rep_dns = tdns;
      sym = sym_ptr[i];
    }
  }
  return {rep_ups, rep_dns, sym};
}

template <class bit_t, class GroupAction>
idx_t ElectronSymmetric<bit_t, GroupAction>::index(bit_t ups, bit_t dns) const {
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

template <class bit_t, class GroupAction>
idx_t ElectronSymmetric<bit_t, GroupAction>::index_switch(bit_t ups,
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

template <class bit_t, class GroupAction>
bool ElectronSymmetric<bit_t, GroupAction>::operator==(
    ElectronSymmetric<bit_t, GroupAction> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) &&
         (charge_conserved_ == rhs.charge_conserved_) &&
         (charge_ == rhs.charge_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}

template <class bit_t, class GroupAction>
bool ElectronSymmetric<bit_t, GroupAction>::operator!=(
    ElectronSymmetric<bit_t, GroupAction> const &rhs) const {
  return !operator==(rhs);
}

template class ElectronSymmetric<uint16, PermutationGroupAction>;
template class ElectronSymmetric<uint32, PermutationGroupAction>;
template class ElectronSymmetric<uint64, PermutationGroupAction>;

} // namespace hydra
