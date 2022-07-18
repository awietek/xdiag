#include "electron_symmetric_indexing.h"

#include <hydra/blocks/utils/block_utils.h>

#include <hydra/indexing/combinations_indexing.h>

#include <hydra/symmetries/fermi_bool_table.h>
#include <hydra/symmetries/representative_list.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra::indexing {

template <class bit_t>
ElectronSymmetricIndexing<bit_t>::ElectronSymmetricIndexing(
    int n_sites, int nup, int ndn, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), n_up_(nup), n_dn_(ndn),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep),
      raw_ups_size_(combinatorics::binomial(n_sites, nup)),
      raw_dns_size_(combinatorics::binomial(n_sites, ndn)),
      lintable_ups_(n_sites, nup), lintable_dns_(n_sites, ndn) {

  using combinatorics::Combinations;

  utils::check_nup_ndn_electron(n_sites, nup, ndn, "Electron");
  utils::check_n_sites(n_sites, permutation_group);

  fermi_bool_ups_table_ = symmetries::fermi_bool_table(
      Combinations<bit_t>(n_sites, nup), group_action_);
  fermi_bool_dns_table_ = symmetries::fermi_bool_table(
      Combinations<bit_t>(n_sites, ndn), group_action_);

  std::tie(reps_up_, idces_up_, syms_up_, sym_limits_up_) =
      symmetries::representatives_indices_symmetries_limits<bit_t>(
          CombinationsIndexing<bit_t>(n_sites, nup), group_action_);

  std::tie(dns_storage_, norms_storage_, dns_limits_, ups_offset_, size_) =
      symmetries::electron_dns_norms_limits_offset_size(
          reps_up_, Combinations<bit_t>(n_sites, ndn), group_action_, irrep_);
}

template class ElectronSymmetricIndexing<uint16_t>;
template class ElectronSymmetricIndexing<uint32_t>;
template class ElectronSymmetricIndexing<uint64_t>;

} // namespace hydra::indexing
