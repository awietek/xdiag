#include "electron_symmetric_v2.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/symmetry_utils.h>

namespace hydra {

template <class bit_t, class GroupAction, class LinTable>
void fill_reps_idces_syms_limits(
    int n_sites, int n_par, GroupAction &&group_action, LinTable &&lintable,
    std::vector<bit_t> &reps, std::vector<idx_t> &idces, std::vector<int> &syms,
    std::vector<std::pair<idx_t, idx_t>> &sym_limits) {

  // Compute all representatives
  idx_t idx = 0;
  for (bit_t state : Combinations(n_sites, n_par)) {
    bit_t rep = utils::representative(state, group_action);
    if (rep == state) {
      reps.push_back(rep);
      idces[idx] = reps.size();
    }
    ++idx;
  }

  // Compute indices of up-representatives and stabilizer symmetries
  idx = 0;
  for (bit_t state : Combinations(n_sites, n_par)) {
    bit_t rep = utils::representative(state, group_action);
    idces[idx] = idces[lintable.index(rep)];

    // Determine the symmetries that yield the up-representative
    idx_t begin = syms.size();
    for (int sym = 0; sym < group_action.n_symmetries(); ++sym) {
      if (group_action.apply(sym, state) == rep)
        syms.push_back(sym);
    }
    idx_t end = syms.size();
    sym_limits[idx] = {begin, end};
    ++idx;
  }
}

template <class bit_t, class GroupAction, class LinTable>
idx_t fill_states_norms(
    int n_sites, int nup, int ndn, GroupAction &&group_action,
    LinTable &&lintable_ups, LinTable &&lintable_dns,
    std::vector<bit_t> const &reps_up, std::vector<idx_t> const &idces_up,
    std::vector<int> const &syms_up,
    std::vector<std::pair<idx_t, idx_t>> const &sym_limits_up,
    std::vector<bool> const &fermi_bool_ups_table,
    std::vector<bool> const &fermi_bool_dns_table, Representation const &irrep,
    std::vector<bit_t> &dns_full, std::vector<double> &norms_dns_full,
    std::map<bit_t, std::vector<bit_t>> &dns_for_up_rep,
    std::map<bit_t, std::vector<double>> &norms_for_up_rep) {

  idx_t raw_ups_size = combinatorics::binomial(n_sites, nup);
  idx_t raw_dns_size = combinatorics::binomial(n_sites, ndn);

  idx_t size = 0;

  // Compute states without stabilizers
  idx_t idx = 0;
  for (bit_t dns : Combinations(n_sites, ndn)) {
    dns_full[idx] = dns;
    norms_dns_full[idx++] = 1.0;
  }

  for (bit_t ups : reps_up) {

    // Get the symmetries that stabilize the ups
    auto [sym_lower, sym_upper] = sym_limits_up[lintable_ups.index(ups)];

    // Only consider ups with non-trivial stabilizer
    if ((sym_upper - sym_lower) > 1) {
      std::vector<bit_t> dn_reps;
      std::vector<double> norms_dn_reps;

      for (bit_t dns : Combinations(n_sites, ndn)) {

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
    } else {
      size += raw_dns_size;
    }
  }
  return size;
}

template <class bit_t, class GroupAction>
ElectronSymmetricV2<bit_t, GroupAction>::ElectronSymmetricV2(
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

  // Compute fermi signs of symmetries acting on all up states
  std::vector<int> fermi_work(n_sites_, 0);
  const int *sym_ptr = group_action_.permutation_array().data();
  for (int sym = 0; sym < n_symmetries_; ++sym) {
    idx_t idx_up = 0;
    for (bit_t ups : Combinations(n_sites, nup)) {
      fermi_bool_ups_table_[sym * raw_ups_size_ + idx_up] =
          utils::fermi_bool_of_permutation(ups, sym_ptr, fermi_work.data());
      ++idx_up;
    }
    sym_ptr += n_sites;
  }

  // Compute fermi signs of symmetries acting on all dn states
  sym_ptr = group_action_.permutation_array().data();
  for (int sym = 0; sym < n_symmetries_; ++sym) {
    idx_t idx_dn = 0;
    for (bit_t dns : Combinations(n_sites, ndn)) {
      fermi_bool_dns_table_[sym * raw_dns_size_ + idx_dn] =
          utils::fermi_bool_of_permutation(dns, sym_ptr, fermi_work.data());
      ++idx_dn;
    }
    sym_ptr += n_sites;
  }

  fill_reps_idces_syms_limits(n_sites, nup, group_action_, lintable_ups_,
                              reps_up_, idces_up_, syms_up_, sym_limits_up_);
  fill_reps_idces_syms_limits(n_sites, ndn, group_action_, lintable_dns_,
                              reps_dn_, idces_dn_, syms_dn_, sym_limits_dn_);

  idx_t size_ups = fill_states_norms(
      n_sites, nup, ndn, group_action_, lintable_ups_, lintable_dns_, reps_up_,
      idces_up_, syms_up_, sym_limits_up_, fermi_bool_ups_table_,
      fermi_bool_dns_table_, irrep_, dns_full_, norms_dns_full_,
      dns_for_up_rep_, norms_for_up_rep_);

  idx_t size_dns = fill_states_norms(
      n_sites, ndn, nup, group_action_, lintable_dns_, lintable_ups_, reps_dn_,
      idces_dn_, syms_dn_, sym_limits_dn_, fermi_bool_dns_table_,
      fermi_bool_ups_table_, irrep_, ups_full_, norms_ups_full_,
      ups_for_dn_rep_, norms_for_dn_rep_);

  assert(size_ups == size_dns);
  size_ = size_ups;
}

template class ElectronSymmetricV2<uint16, PermutationGroupAction>;
template class ElectronSymmetricV2<uint32, PermutationGroupAction>;
template class ElectronSymmetricV2<uint64, PermutationGroupAction>;

template class ElectronSymmetricV2<uint16, PermutationGroupLookup<uint16>>;
template class ElectronSymmetricV2<uint32, PermutationGroupLookup<uint32>>;
template class ElectronSymmetricV2<uint64, PermutationGroupLookup<uint64>>;

} // namespace hydra
