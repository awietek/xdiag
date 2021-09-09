#include "indexing_symmetric.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/models/utils/model_utils.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/symmetry_utils.h>

namespace hydra::indexing {

  template <class bit_t, class GroupAction>
IndexingSymmetric<bit_t, GroupAction>::IndexingSymmetric(
    int n_sites, int n_up, PermutationGroup permutation_group,
    Representation irrep)
    : lin_table_(n_sites, n_up),
      index_of_raw_state_(combinatorics::binomial(n_sites, n_up),
                          invalid_index),
      norm_of_raw_state_(combinatorics::binomial(n_sites, n_up), 0.) {
  utils::check_nup_spinhalf(n_sites, n_up, "IndexingSymmetric");

  // if not all symmetries are allowed by irrep, choose a subgroup
  if (irrep.allowed_symmetries().size() > 0) {
    permutation_group = permutation_group.subgroup(irrep.allowed_symmetries());
  }
  auto group_action = GroupAction(permutation_group);

  // Go through non symmetrized states and register representatives
  idx_t idx = 0;
  idx_t n_representatives = 0;
  
  for (bit_t state : Combinations(n_sites, n_up)) {

    bit_t rep = utils::representative(state, group_action);
    // register state if it's a representative
    if (rep == state) {

      double norm = utils::compute_norm(rep, group_action, irrep);

      // if (norm > 1e-6) { // tolerance big as 1e-6 since root is taken
        idx_t idx = lin_table_.index(rep);
        states_.push_back(rep);
        norm_of_raw_state_[idx] = norm;
        index_of_raw_state_[idx] = n_representatives++;
      // }
    }
    ++idx;
  }

  // Go through non symmetrized states and fill indices for all states
  idx = 0;
  for (bit_t state : Combinations(n_sites, n_up)) {

    // Compute the representative of state and the corresponding symmetry
    auto [rep, rep_sym] = utils::representative_sym(state, group_action);

    if (rep != state) {
      complex norm = norm_of_raw_state_[lin_table_.index(rep)];
      // if (std::abs(norm) > 1e-6) {
        norm_of_raw_state_[idx] = irrep.character(rep_sym) * norm;
        index_of_raw_state_[idx] = index_of_raw_state_[lin_table_.index(rep)];
      // }
    }
    ++idx;
  }

  size_ = (idx_t)states_.size();
}

template class IndexingSymmetric<uint16, PermutationGroupAction>;
template class IndexingSymmetric<uint32, PermutationGroupAction>;
template class IndexingSymmetric<uint64, PermutationGroupAction>;

template class IndexingSymmetric<uint16, PermutationGroupLookup<uint16>>;
template class IndexingSymmetric<uint32, PermutationGroupLookup<uint32>>;
template class IndexingSymmetric<uint64, PermutationGroupLookup<uint64>>;


} // namespace hydra::indexing
