#include "electron_symmetric.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra {

template <class bit_t, class GroupAction, class LinTable>
idx_t fill_states_norms_electron(
    int n_sites, int nup, int ndn, GroupAction &&group_action,
    LinTable &&lintable_ups, LinTable &&lintable_dns,
    std::vector<bit_t> const &reps_up, std::vector<idx_t> const &idces_up,
    std::vector<int> const &syms_up,
    std::vector<std::pair<gsl::span<int const>::size_type,
                          gsl::span<int const>::size_type>> const
        &sym_limits_up,
    std::vector<bool> const &fermi_bool_ups_table,
    std::vector<bool> const &fermi_bool_dns_table, Representation const &irrep,
    std::vector<idx_t> &up_offsets, std::vector<bit_t> &dns_full,
    std::vector<double> &norms_dns_full,
    std::map<bit_t, std::vector<bit_t>> &dns_for_up_rep,
    std::map<bit_t, std::vector<double>> &norms_for_up_rep) {

  using combinatorics::Combinations;

  // create full dns states, used when ups have trivial stabilizer
  idx_t idx = 0;
  for (bit_t dns : Combinations<bit_t>(n_sites, ndn)) {
    dns_full[idx] = dns;
    norms_dns_full[idx++] = 1.0;
  }

  // create the dns states, when ups have non-trivial stabilizer
  idx_t size = 0;
  idx_t idx_up = 0;
  up_offsets.resize(reps_up.size());
  for (bit_t ups : reps_up) {

    up_offsets[idx_up] = size;

    // Get the symmetries that stabilize the ups
    auto [start, length] = sym_limits_up[lintable_ups.index(ups)];  // HOTFIX
    idx_t sym_lower = start;
    idx_t sym_upper = start + length;
    // auto [sym_lower, sym_upper] = sym_limits_up[lintable_ups.index(ups)];

    // ups have trivial stabilizer, we don't store dns, but use dns_full
    if ((sym_upper - sym_lower) == 1) {
      size += combinatorics::binomial(n_sites, ndn);
    }

    // ups have non-trivial stabilizer, we store the dns configurations
    else {

      std::vector<bit_t> dn_reps;
      std::vector<double> norms_dn_reps;
      std::vector<int> syms_stable(syms_up.begin() + sym_lower,
                                   syms_up.begin() + sym_upper);

      for (bit_t dns : Combinations<bit_t>(n_sites, ndn)) {

        bit_t dn_rep =
            symmetries::representative_subset(dns, group_action, syms_stable);

        if (dns == dn_rep) { // if "dns" is representative

          double norm = symmetries::norm_electron_subset(ups, dns, group_action,
                                                         irrep, syms_stable);

          if (norm > 1e-6) { // only keep dns with non-zero norm
            dn_reps.push_back(dn_rep);
            norms_dn_reps.push_back(norm);
          }
        }
        dns_for_up_rep[ups] = dn_reps;
        norms_for_up_rep[ups] = norms_dn_reps;
      }
      size += dn_reps.size();
    }
    ++idx_up;
  }
  return size;
}

template <class bit_t, class GroupAction>
ElectronSymmetric<bit_t, GroupAction>::ElectronSymmetric(
    int n_sites, int nup, int ndn, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), charge_conserved_(true), charge_(nup + ndn),
      sz_conserved_(true), sz_(nup - ndn), n_up_(nup), n_dn_(ndn),
      permutation_group_(
          (irrep.allowed_symmetries().size() > 0)
              ? permutation_group.subgroup(irrep.allowed_symmetries())
              : permutation_group),
      group_action_(permutation_group_), irrep_(irrep),
      n_symmetries_(group_action_.n_symmetries()), lintable_ups_(n_sites, nup),
      lintable_dns_(n_sites, ndn),
      raw_ups_size_(combinatorics::binomial(n_sites, nup)),
      raw_dns_size_(combinatorics::binomial(n_sites, ndn)),
      fermi_bool_ups_table_(
          symmetries::fermi_bool_table<bit_t>(nup, group_action_)),
      fermi_bool_dns_table_(
          symmetries::fermi_bool_table<bit_t>(ndn, group_action_)),
      dns_full_(raw_dns_size_), norms_dns_full_(raw_dns_size_),
      ups_full_(raw_ups_size_), norms_ups_full_(raw_ups_size_) {

  utils::check_nup_ndn_electron(n_sites, nup, ndn, "ElectronSymmetric");

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          nup, group_action_, lintable_ups_);

  std::tie(reps_dn_, idces_dn_, syms_dn_, sym_limits_dn_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          ndn, group_action_, lintable_dns_);

  idx_t size_ups = fill_states_norms_electron(
      n_sites, nup, ndn, group_action_, lintable_ups_, lintable_dns_, reps_up_,
      idces_up_, syms_up_, sym_limits_up_, fermi_bool_ups_table_,
      fermi_bool_dns_table_, irrep_, up_offsets_, dns_full_, norms_dns_full_,
      dns_for_up_rep_, norms_for_up_rep_);

  idx_t size_dns = fill_states_norms_electron(
      n_sites, ndn, nup, group_action_, lintable_dns_, lintable_ups_, reps_dn_,
      idces_dn_, syms_dn_, sym_limits_dn_, fermi_bool_dns_table_,
      fermi_bool_ups_table_, irrep_, dn_offsets_, ups_full_, norms_ups_full_,
      ups_for_dn_rep_, norms_for_dn_rep_);

  // lila::Log("size_ups: {}, size_dns: {}", size_ups, size_dns);
  assert(size_ups == size_dns);
  size_ = size_ups;
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

template class ElectronSymmetric<uint16_t, PermutationGroupAction>;
template class ElectronSymmetric<uint32_t, PermutationGroupAction>;
template class ElectronSymmetric<uint64_t, PermutationGroupAction>;

template class ElectronSymmetric<uint16_t, PermutationGroupLookup<uint16_t>>;
template class ElectronSymmetric<uint32_t, PermutationGroupLookup<uint32_t>>;
template class ElectronSymmetric<uint64_t, PermutationGroupLookup<uint64_t>>;

} // namespace hydra
