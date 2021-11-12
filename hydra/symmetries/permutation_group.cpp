#include "permutation_group.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include <hydra/symmetries/symmetry_operations.h>
#include <lila/utils/logger.h>

namespace hydra {

PermutationGroup::PermutationGroup(
    std::vector<std::vector<int>> const &symmetries)
    : n_sites_(symmetries[0].size()), n_symmetries_(symmetries.size()),
      permutation_array_(n_sites_ * n_symmetries_) {
  for (int i = 0; i < n_symmetries_; ++i) {

    // Check whether lattice symmetries are well-formed
    if ((int)symmetries[i].size() != n_sites_)
      lila::Log.err("Error constructing PermutationGroup: "
                    "there's a symmetry not of length n_sites");

    if (!symmetries::is_valid_permutation(n_sites_, symmetries[i].data()))
      lila::Log.err("Error constructing PermutationGroup: "
                    "there's a symmetry which is not a proper permutation");

    std::copy(symmetries[i].begin(), symmetries[i].end(),
              permutation_array_.begin() + i * n_sites_);
  }
}

PermutationGroup::PermutationGroup(int n_sites, int n_symmetries,
                                   std::vector<int> const &permutation_array)
    : n_sites_(n_sites), n_symmetries_(n_symmetries),
      permutation_array_(permutation_array) {

  if ((int)permutation_array.size() != n_sites * n_symmetries)
    lila::Log.err("Error constructing PermutationGroup: "
                  "there's a symmetry not of length n_sites");
  for (int i = 0; i < n_symmetries; ++i) {
    if (!symmetries::is_valid_permutation(n_sites_,
                                     permutation_array_.data() + i * n_sites_))
      lila::Log.err("Error constructing PermutationGroup: "
                    "there's a symmetry which is not a proper permutation");
  }
}

PermutationGroup
PermutationGroup::subgroup(std::vector<int> const &symmetry_numbers) const {
  std::vector<int> subgroup_permutation_array;

  for (int n_sym : symmetry_numbers) {
    if ((0 > n_sym) || (n_sym >= n_symmetries_))
      lila::Log.err("Error building subgroup of PermutationGroup: "
                    "invalid symmetry index");
    subgroup_permutation_array.insert(
        subgroup_permutation_array.end(),
        permutation_array_.begin() + n_sym * n_sites_,
        permutation_array_.begin() + (n_sym + 1) * n_sites_);
  }

  return PermutationGroup(n_sites_, symmetry_numbers.size(),
                          subgroup_permutation_array);
}

bool PermutationGroup::operator==(PermutationGroup const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_symmetries_ == rhs.n_symmetries_) &&
         (permutation_array_ == rhs.permutation_array_);
}

bool PermutationGroup::operator!=(PermutationGroup const &rhs) const {
  return !operator==(rhs);
}

} // namespace hydra
