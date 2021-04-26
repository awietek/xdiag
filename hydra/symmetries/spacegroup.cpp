#include <cassert>
#include <fstream>
#include <iostream>

#include "spacegroup.h"

namespace hydra {

template <class bit_t, class SpaceGroupOperator>
SpaceGroup<bit_t, SpaceGroupOperator>::SpaceGroup(
    std::vector<std::vector<int>> const &symmetries)
    : n_sites_(symmetries[0].size()), n_symmetries_(symmetries.size()),
      symmetries_(symmetries), permutation_array_(n_sites_ * n_symmetries_) {
  for (int i = 0; i < (int)symmetries.size(); ++i) {

    // Check whether lattice symmetries are well-formed
    if ((int)symmetries[i].size() != n_sites_)
      HydraLog.err("Error constructing SpaceGroup: "
                   "there's a symmetry not of length n_sites");

    if (!symmetries::is_valid_permutation(n_sites_, symmetries[i].data()))
      HydraLog.err("Error constructing SpaceGroup: "
                   "there's a symmetry which is not a proper permutation");

    std::copy(symmetries[i].begin(), symmetries[i].end(),
              permutation_array_.begin() + i * n_sites_);
  }
  spacegroup_operator_ = SpaceGroupOperator(n_sites_, permutation_array_);
}

template <class bit_t, class SpaceGroupOperator>
SpaceGroup<bit_t, SpaceGroupOperator>
SpaceGroup<bit_t, SpaceGroupOperator>::subgroup(
    std::vector<int> const &symmetry_numbers) const {
  std::vector<std::vector<int>> subgroup_symmetries;
  for (int n_sym : symmetry_numbers) {
    if ((0 > n_sym) && (n_sym >= (int)symmetries_.size()))
      HydraLog.err("Error building subgroup of SpaceGroup: "
                   "invalid symmetry index");
    subgroup_symmetries.push_back(symmetries_[n_sym]);
  }
  return SpaceGroup<bit_t, SpaceGroupOperator>(subgroup_symmetries);
}

template <class bit_t, class SpaceGroupOperator>
bool SpaceGroup<bit_t, SpaceGroupOperator>::operator==(
    SpaceGroup const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_symmetries_ == rhs.n_symmetries_) &&
         (spacegroup_operator_ == rhs.spacegroup_operator_);
}
template <class bit_t, class SpaceGroupOperator>
bool SpaceGroup<bit_t, SpaceGroupOperator>::operator!=(
    SpaceGroup const &rhs) const {
  return !operator==(rhs);
}

template class SpaceGroup<uint16, SpaceGroupOperator<uint16>>;
template class SpaceGroup<uint32, SpaceGroupOperator<uint32>>;
template class SpaceGroup<uint64, SpaceGroupOperator<uint64>>;

} // namespace hydra
