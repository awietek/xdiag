#include "spinhalf_symmetric_indexing.h"

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/symmetries/representative_list.h>

namespace hydra::indexing {

template <class bit_t>
SpinhalfSymmetricIndexing<bit_t>::SpinhalfSymmetricIndexing(
    int n_sites, int n_up, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), n_up_(n_up),
      group_action_((irrep.allowed_symmetries().size() > 0)
                        ? permutation_group.subgroup(irrep.allowed_symmetries())
                        : permutation_group),
      irrep_(irrep), lin_table_(n_sites, n_up) {

  utils::check_nup_spinhalf(n_sites, n_up, "SpinhalfSymmetricIndexing");
  utils::check_n_sites(n_sites, permutation_group);

  std::tie(reps_, index_for_rep_, syms_, sym_limits_for_rep_, norms_) =
      symmetries::representatives_indices_symmetries_limits_norms<bit_t>(
          group_action_, lin_table_, irrep);
  size_ = (idx_t)reps_.size();
}

template class SpinhalfSymmetricIndexing<uint16_t>;
template class SpinhalfSymmetricIndexing<uint32_t>;
template class SpinhalfSymmetricIndexing<uint64_t>;

} // namespace hydra::indexing