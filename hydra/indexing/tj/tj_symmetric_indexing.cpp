#include "tj_symmetric_indexing.h"

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/indexing/combinations_indexing.h>

#include <hydra/combinatorics/combinations.h>
#include <hydra/symmetries/operations/symmetry_operations.h>

namespace hydra::indexing {

template <typename bit_t>
tJSymmetricIndexing<bit_t>::tJSymmetricIndexing(int n_sites, int nup, int ndn,
                                                PermutationGroup group,
                                                Representation irrep)
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn),
      group_action_(allowed_subgroup(group, irrep)), irrep_(irrep),
      raw_ups_size_(combinatorics::binomial(n_sites, nup)),
      raw_dns_size_(combinatorics::binomial(n_sites, ndn)),
      raw_dnsc_size_(combinatorics::binomial(n_sites - nup, ndn)),
      lintable_ups_(n_sites, nup), lintable_dns_(n_sites, ndn),
      lintable_dnsc_(n_sites - nup, ndn), fermi_table_ups_(n_sites, nup, group),
      fermi_table_dns_(n_sites, ndn, group) {

  using combinatorics::Combinations;
  utils::check_n_sites(n_sites, group);
  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          lintable_ups_, group_action_);

  // if ups have trivial stabilizer, dns (compressed) are stored in front
  for (bit_t dns : Combinations<bit_t>(n_sites - nup, ndn)) {
    dns_storage_.push_back(dns);
    norms_storage_.push_back(1.0);
  }

  ups_offset_.resize(reps_up_.size());
  dns_limits_.resize(reps_up_.size());
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  size_ = 0;
  idx_t idx_up = 0;

  for (bit_t ups : reps_up_) {
    ups_offset_[idx_up] = size_;
    auto syms = syms_ups(ups);

    // ups have trivial stabilizer -> dns stored in beginning
    if (syms.size() == 1) {
      span_size_t start = 0;
      span_size_t length = raw_dnsc_size_;
      dns_limits_[idx_up] = {start, length};
      size_ += length;
      // ups have non-trivial stabilizer, we store the dns configurations
    } else {
      bit_t not_ups = (~ups) & sitesmask;

      span_size_t start = dns_storage_.size();
      for (bit_t dnsc : Combinations<bit_t>(n_sites - nup, ndn)) {
        bit_t dns = bitops::deposit(dnsc, not_ups);
        bit_t dns_rep =
            symmetries::representative_subset(dns, group_action_, syms);
        if (dns == dns_rep) {
          double norm = symmetries::norm_electron_subset(
              ups, dns, group_action_, irrep, syms);
          if (norm > 1e-6) { // only keep dns with non-zero norm
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

template class tJSymmetricIndexing<uint16_t>;
template class tJSymmetricIndexing<uint32_t>;
template class tJSymmetricIndexing<uint64_t>;

} // namespace hydra::indexing
