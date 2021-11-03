#include "electron_symmetric.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/symmetry_utils.h>

namespace hydra {

template <class bit_t, class GroupAction, class LinTable>
idx_t fill_states_norms_electron(
    int n_sites, int nup, int ndn, GroupAction &&group_action,
    LinTable &&lintable_ups, LinTable &&lintable_dns,
    std::vector<bit_t> const &reps_up, std::vector<idx_t> const &idces_up,
    std::vector<int> const &syms_up,
    std::vector<std::pair<idx_t, idx_t>> const &sym_limits_up,
    std::vector<bool> const &fermi_bool_ups_table,
    std::vector<bool> const &fermi_bool_dns_table, Representation const &irrep,
    std::vector<idx_t> &up_offsets, std::vector<bit_t> &dns_full,
    std::vector<double> &norms_dns_full,
    std::map<bit_t, std::vector<bit_t>> &dns_for_up_rep,
    std::map<bit_t, std::vector<double>> &norms_for_up_rep) {

  using combinatorics::Combinations;

  idx_t raw_ups_size = combinatorics::binomial(n_sites, nup);
  idx_t raw_dns_size = combinatorics::binomial(n_sites, ndn);

  // Compute states without stabilizers
  idx_t idx = 0;
  for (bit_t dns : Combinations<bit_t>(n_sites, ndn)) {
    dns_full[idx] = dns;
    norms_dns_full[idx++] = 1.0;
  }

  idx_t size = 0;
  idx_t idx_up = 0;
  up_offsets.resize(reps_up.size());
  for (bit_t ups : reps_up) {

    up_offsets[idx_up] = size;

    // Get the symmetries that stabilize the ups
    auto [sym_lower, sym_upper] = sym_limits_up[lintable_ups.index(ups)];

    if ((sym_upper - sym_lower) == 1) {
      size += raw_dns_size;
    } else {
      std::vector<bit_t> dn_reps;
      std::vector<double> norms_dn_reps;

      for (bit_t dns : Combinations<bit_t>(n_sites, ndn)) {

        // Determine dn representative
        bit_t dn_rep = dns;
        for (idx_t sym_idx = sym_lower; sym_idx < sym_upper; ++sym_idx) {
          int sym = syms_up[sym_idx];
          bit_t tdns = group_action.apply(sym, dns);
          if (tdns < dn_rep)
            dn_rep = tdns;
        }

        // if "dns" is representative
        if (dns == dn_rep) {

          // Compute its norm ...
          complex amplitude = 0.0;
          for (idx_t sym_idx = sym_lower; sym_idx < sym_upper; ++sym_idx) {
            int sym = syms_up[sym_idx];
            assert(group_action.apply(sym, ups) == ups);
            if (group_action.apply(sym, dn_rep) == dn_rep) {
              bool fermi_up = fermi_bool_ups_table[sym * raw_ups_size +
                                                   lintable_ups.index(ups)];
              bool fermi_dn = fermi_bool_dns_table[sym * raw_dns_size +
                                                   lintable_dns.index(dns)];
              if (fermi_up == fermi_dn) {
                amplitude += irrep.character(sym);
              } else {
                amplitude -= irrep.character(sym);
              }
            }
          }
          double norm = std::sqrt(std::abs(amplitude));

          // ... and keep state if norm non-zero
          if (norm > 1e-6) {
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
      fermi_bool_ups_table_(raw_ups_size_ * n_symmetries_),
      fermi_bool_dns_table_(raw_dns_size_ * n_symmetries_),
      idces_up_(raw_ups_size_), sym_limits_up_(raw_ups_size_, {0, 0}),
      idces_dn_(raw_dns_size_), sym_limits_dn_(raw_dns_size_, {0, 0}),
      dns_full_(raw_dns_size_), norms_dns_full_(raw_dns_size_),
      ups_full_(raw_ups_size_), norms_ups_full_(raw_ups_size_) {

  using combinatorics::Combinations;

  utils::fill_fermi_bool_table<bit_t>(group_action_, nup, fermi_bool_ups_table_);
  utils::fill_fermi_bool_table<bit_t>(group_action_, ndn, fermi_bool_dns_table_);

  utils::fill_reps_idces_syms_limits(n_sites, nup, group_action_, lintable_ups_,
                                     reps_up_, idces_up_, syms_up_,
                                     sym_limits_up_);
  utils::fill_reps_idces_syms_limits(n_sites, ndn, group_action_, lintable_dns_,
                                     reps_dn_, idces_dn_, syms_dn_,
                                     sym_limits_dn_);

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
