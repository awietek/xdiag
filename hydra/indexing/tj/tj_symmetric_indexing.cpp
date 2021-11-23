#include "tj_symmetric_indexing.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra::indexing {

template <typename bit_t, class GroupAction>
static idx_t fill_states_norms_tj(
    int n_sites, int nup, int ndn, GroupAction &&group_action,
    LinTable<bit_t> const &lintable_ups, LinTable<bit_t> const &lintable_dns,
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

  // create full dns states, used when ups have trivial stabilizer
  // for the tJ block, dns_full are compressed on not_ups sites
  idx_t idx = 0;
  for (bit_t dns : Combinations<bit_t>(n_sites - nup, ndn)) {
    dns_full[idx] = dns;
    norms_dns_full[idx++] = 1.0;
  }

  // create the dns states, when ups have non-trivial stabilizer_symmetries
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;
  idx_t size = 0;
  idx_t idx_up = 0;
  up_offsets.resize(reps_up.size());
  for (bit_t ups : reps_up) {
    up_offsets[idx_up] = size;

    // Get the symmetries that stabilize the ups
    auto [sym_lower, sym_upper] = sym_limits_up[lintable_ups.index(ups)];

    bit_t not_ups = (~ups) & sitesmask;

    // ups have trivial stabilizer, we don't store dns, but use dns_full
    if ((sym_upper - sym_lower) == 1) {
      size += combinatorics::binomial(n_sites - nup, ndn);
    }

    // ups have non-trivial stabilizer, we store the dns configurations
    else {
      std::vector<bit_t> dn_reps;
      std::vector<double> norms_dn_reps;
      std::vector<int> syms_stable(syms_up.begin() + sym_lower,
                                   syms_up.begin() + sym_upper);

      for (bit_t dns_compressed : Combinations<bit_t>(n_sites - nup, ndn)) {
        bit_t dns = bitops::deposit(dns_compressed, not_ups);

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

template <class bit_t>
tJSymmetricIndexing<bit_t>::tJSymmetricIndexing(
    int n_sites, int nup, int ndn, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn),
      group_action_(permutation_group), irrep_(irrep),
      raw_ups_size_(combinatorics::binomial(n_sites, nup)),
      raw_dns_size_(combinatorics::binomial(n_sites, ndn)),
      raw_upsc_size_(combinatorics::binomial(n_sites - ndn, nup)),
      raw_dnsc_size_(combinatorics::binomial(n_sites - nup, ndn)),
      lintable_ups_(n_sites, nup), lintable_dns_(n_sites, ndn),
      lintable_upsc_(n_sites - ndn, nup), lintable_dnsc_(n_sites - nup, ndn),
      dns_full_(raw_dnsc_size_), norms_dns_full_(raw_dnsc_size_)
      // ups_full_(raw_upsc_size_), norms_ups_full_(raw_upsc_size_)
{

  PermutationGroupLookup<bit_t> group_action(permutation_group);
  fermi_bool_ups_table_ =
      symmetries::fermi_bool_table<bit_t>(nup, group_action);
  fermi_bool_dns_table_ =
      symmetries::fermi_bool_table<bit_t>(ndn, group_action);

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          nup, group_action, lintable_ups_);
  // std::tie(reps_dn_, idces_dn_, syms_dn_, sym_limits_dn_) =
  //     symmetries::representatives_indices_symmetries_limits<bit_t>(
  //         ndn, group_action, lintable_dns_);

  idx_t size_ups = fill_states_norms_tj(
      n_sites, nup, ndn, group_action, lintable_ups_, lintable_dns_, reps_up_,
      idces_up_, syms_up_, sym_limits_up_, fermi_bool_ups_table_,
      fermi_bool_dns_table_, irrep, up_offsets_, dns_full_, norms_dns_full_,
      dns_for_up_rep_, norms_for_up_rep_);

  // idx_t size_dns = fill_states_norms_tj(
  //     n_sites, ndn, nup, group_action, lintable_dns_, lintable_ups_, reps_dn_,
  //     idces_dn_, syms_dn_, sym_limits_dn_, fermi_bool_dns_table_,
  //     fermi_bool_ups_table_, irrep, dn_offsets_, ups_full_, norms_ups_full_,
  //     ups_for_dn_rep_, norms_for_dn_rep_);

  size_ = size_ups;
}

template class tJSymmetricIndexing<uint16_t>;
template class tJSymmetricIndexing<uint32_t>;
template class tJSymmetricIndexing<uint64_t>;

} // namespace hydra::indexing
