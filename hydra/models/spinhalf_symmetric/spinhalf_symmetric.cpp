#include "spinhalf_symmetric.h"

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>

#include <hydra/models/utils/model_utils.h>
#include <hydra/models/utils/symmetrized_norm.h>

namespace hydra {

template <class bit_t, class GroupAction>
SpinhalfSymmetric<bit_t, GroupAction>::SpinhalfSymmetric(
    int n_sites, int n_up, PermutationGroup permutation_group,
    Representation irrep)
    : n_sites_(n_sites), sz_conserved_(true), n_up_(n_up),
      n_dn_(n_sites - n_up), sz_(n_up_ - n_dn_),
      permutation_group_(permutation_group), irrep_(irrep),
      indexing_(n_sites, n_up, permutation_group, irrep) {
  utils::check_nup_spinhalf(n_sites, n_up, "SpinhalfSymmetric");

  // if not all symmetries are allowed by irrep, choose a subgroup
  if (irrep.allowed_symmetries().size() > 0) {
    permutation_group_ = permutation_group.subgroup(irrep.allowed_symmetries());
  }
}

template <class bit_t, class GroupAction>
bool SpinhalfSymmetric<bit_t, GroupAction>::operator==(
    SpinhalfSymmetric<bit_t, GroupAction> const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (sz_conserved_ == rhs.sz_conserved_) &&
         (sz_ == rhs.sz_) && (n_up_ == rhs.n_up_) && (n_dn_ == rhs.n_dn_) &&
         (permutation_group_ == rhs.permutation_group_) &&
         (irrep_ == rhs.irrep_);
}

template <class bit_t, class GroupAction>
bool SpinhalfSymmetric<bit_t, GroupAction>::operator!=(
    SpinhalfSymmetric<bit_t, GroupAction> const &rhs) const {
  return !operator==(rhs);
}

template class SpinhalfSymmetric<uint16, PermutationGroupAction>;
template class SpinhalfSymmetric<uint32, PermutationGroupAction>;
template class SpinhalfSymmetric<uint64, PermutationGroupAction>;

} // namespace hydra
