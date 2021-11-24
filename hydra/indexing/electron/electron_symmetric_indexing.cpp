#include "electron_symmetric_indexing.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra::indexing {

template <class bit_t>
ElectronSymmetricIndexing<bit_t>::ElectronSymmetricIndexing(
    int n_sites, int nup, int ndn, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn),
      group_action_(permutation_group), irrep_(irrep),
      raw_ups_size_(combinatorics::binomial(n_sites, nup)),
      raw_dns_size_(combinatorics::binomial(n_sites, ndn)),
      lintable_ups_(n_sites, nup), lintable_dns_(n_sites, ndn) {

  using combinatorics::Combinations;

  PermutationGroupLookup<bit_t> group_action(permutation_group);
  fermi_bool_ups_table_ =
      symmetries::fermi_bool_table<bit_t>(nup, group_action);
  fermi_bool_dns_table_ =
      symmetries::fermi_bool_table<bit_t>(ndn, group_action);

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          nup, group_action, lintable_ups_);

  // if ups have trivial stabilizer, dns  are stored in front
  for (bit_t dns : Combinations<bit_t>(n_sites, ndn)) {
    dns_storage_.push_back(dns);
    norms_storage_.push_back(1.0);
  }

  ups_offset_.resize(reps_up_.size());
  dns_limits_.resize(reps_up_.size());

  size_ = 0;
  idx_t idx_up = 0;

  for (bit_t ups : reps_up_) {
    ups_offset_[idx_up] = size_;
    auto syms = syms_ups(ups);

    // ups have trivial stabilizer -> dns stored in beginning
    if (syms.size() == 1) {
      span_size_t start = 0;
      span_size_t length = raw_dns_size_;
      dns_limits_[idx_up] = {start, length};
      size_ += length;
      // ups have non-trivial stabilizer, we store the dns configurations
    } else {
      span_size_t start = dns_storage_.size();
      for (bit_t dns : Combinations<bit_t>(n_sites, ndn)) {
        bit_t dns_rep =
            symmetries::representative_subset(dns, group_action, syms);
        if (dns == dns_rep) {
          double norm = symmetries::norm_electron_subset(ups, dns, group_action,
                                                         irrep, syms);
          if (norm > 1e-6) {
            dns_storage_.push_back(dns_rep);
            norms_storage_.push_back(norm);
          }
        }
      }
      span_size_t end = dns_storage_.size();
      span_size_t length = end - start;
      dns_limits_[idx_up] = {start, length};
      size_ += length;
    }
    ++idx_up;
  }
}

template class ElectronSymmetricIndexing<uint16_t>;
template class ElectronSymmetricIndexing<uint32_t>;
template class ElectronSymmetricIndexing<uint64_t>;

} // namespace hydra::indexing
