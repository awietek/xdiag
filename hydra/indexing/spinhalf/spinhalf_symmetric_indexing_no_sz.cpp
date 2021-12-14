#include "spinhalf_symmetric_indexing_no_sz.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/symmetries/representative_list.h>

namespace hydra::indexing {

template <class bit_t>
SpinhalfSymmetricIndexingNoSz<bit_t>::SpinhalfSymmetricIndexingNoSz(
    int n_sites, PermutationGroup permutation_group, Representation irrep)
    : n_sites_(n_sites),
      group_action_(allowed_subgroup(permutation_group, irrep)), irrep_(irrep),
      subsets_indexing_(n_sites) {

  assert(n_sites >= 0);
  utils::check_n_sites(n_sites, permutation_group);

  std::tie(reps_, index_for_rep_, syms_, sym_limits_for_rep_, norms_) =
      symmetries::representatives_indices_symmetries_limits_norms<bit_t>(
          subsets_indexing_, group_action_, irrep);

  size_ = (idx_t)reps_.size();
}

template class SpinhalfSymmetricIndexingNoSz<uint16_t>;
template class SpinhalfSymmetricIndexingNoSz<uint32_t>;
template class SpinhalfSymmetricIndexingNoSz<uint64_t>;

} // namespace hydra::indexing
