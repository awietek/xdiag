#include "indexing_symmetric_fermionic.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/permutation_group_lookup.h>
#include <hydra/symmetries/symmetry_utils.h>

namespace hydra::indexing {

template <class bit_t, class GroupAction>
IndexingSymmetricFermionic<bit_t, GroupAction>::IndexingSymmetricFermionic(
    int n_sites, int n_fermions, PermutationGroup permutation_group,
    Representation irrep)
    : lin_table_(n_sites, n_fermions),
      index_of_raw_state_(combinatorics::binomial(n_sites, n_fermions),
                          invalid_index),
      norm_of_raw_state_(combinatorics::binomial(n_sites, n_fermions), 0.) {
  using combinatorics::Combinations;

  utils::check_nup_spinhalf(n_sites, n_fermions, "IndexingSymmetricFermionic");

  // if not all symmetries are allowed by irrep, choose a subgroup
  if (irrep.allowed_symmetries().size() > 0) {
    permutation_group = permutation_group.subgroup(irrep.allowed_symmetries());
  }
  auto group_action = GroupAction(permutation_group);

  // Go through non symmetrized states and register representatives
  idx_t idx = 0;
  idx_t n_representatives = 0;
  std::vector<int> fermi_work(n_sites, 0);

  for (bit_t state : Combinations(n_sites, n_fermions)) {

    bit_t rep = utils::representative(state, group_action);
    // register state if it's a representative
    if (rep == state) {

      double norm = utils::compute_norm_fermionic(rep, group_action, irrep,
                                                  fermi_work.data());

      if (norm > 1e-6) { // tolerance big as 1e-6 since root is taken
        idx_t idx = lin_table_.index(rep);
        states_.push_back(rep);
        norm_of_raw_state_[idx] = norm;
        index_of_raw_state_[idx] = n_representatives++;
      }
    }
    ++idx;
  }

  // Go through non symmetrized states and fill indices for all states
  idx = 0;
  const int *sym_ptr = permutation_group.permutation_array().data();
  for (bit_t state : Combinations(n_sites, n_fermions)) {

    // Compute the representative of state and the corresponding symmetry
    auto [rep, rep_sym] = utils::representative_sym(state, group_action);

    if (rep != state) {
      complex norm = norm_of_raw_state_[lin_table_.index(rep)];
      if (std::abs(norm) > 1e-6) {
        norm_of_raw_state_[idx] =
            utils::fermi_sign_of_permutation(state, sym_ptr + rep_sym * n_sites,
                                             fermi_work.data()) *
            irrep.character(rep_sym) * norm;
        index_of_raw_state_[idx] = index_of_raw_state_[lin_table_.index(rep)];
      }
    }
    ++idx;
  }

  size_ = (idx_t)states_.size();
}

template class IndexingSymmetricFermionic<uint16, PermutationGroupAction>;
template class IndexingSymmetricFermionic<uint32, PermutationGroupAction>;
template class IndexingSymmetricFermionic<uint64, PermutationGroupAction>;

template class IndexingSymmetricFermionic<uint16,
                                          PermutationGroupLookup<uint16>>;
template class IndexingSymmetricFermionic<uint32,
                                          PermutationGroupLookup<uint32>>;
template class IndexingSymmetricFermionic<uint64,
                                          PermutationGroupLookup<uint64>>;

} // namespace hydra::indexing
