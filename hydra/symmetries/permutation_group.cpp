#include "permutation_group.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include <hydra/symmetries/operations/symmetry_operations.h>
#include <hydra/utils/logger.h>

namespace hydra {

PermutationGroup::PermutationGroup(std::vector<Permutation> const &permutations)
    : n_sites_(permutations.size() > 0 ? permutations[0].n_sites() : 0),
      n_symmetries_(permutations.size()), permutations_(permutations),
      inverse_(n_symmetries_) {

  // Check whether all permutations have same number of sites
  for (auto p : permutations) {
    if (p.n_sites() != n_sites_) {
      Log.err("Error constructing PermutationGroup: not all Permutations have "
              "the same number of sites");
    }
  }

  // Check whether identity is contained
  if (n_sites_ > 0) {
    auto id = identity_permutation(n_sites_);
    if (std::find(permutations.begin(), permutations.end(), id) ==
        permutations.end()) {
      Log.err("Error constructing PermutationGroup: no identity element found");
    }
  }
  
  // Check multiplication is complete
  for (auto p1 : permutations) {
    for (auto p2 : permutations) {
      auto p = p1 * p2;
      if (std::find(permutations.begin(), permutations.end(), p) ==
          permutations.end()) {
        Log.err("Error constructing PermutationGroup: group multiplication not "
                "closed");
      }
    }
  }

  // Check if inverse exists
  int idx = 0;
  for (auto p : permutations) {
    auto pinv = hydra::inverse(p);
    auto it = std::find(permutations.begin(), permutations.end(), pinv);
    if (it == permutations.end()) {
      Log.err("Error constructing PermutationGroup: inverse element not found");
    } else {
      int idx_inv = std::distance(permutations.begin(), it);
      inverse_[idx] = idx_inv;
    }
    idx++;
  }
}

PermutationGroup
PermutationGroup::subgroup(std::vector<int> const &symmetry_numbers) const {
  std::vector<Permutation> subgroup_permutations;

  for (int n_sym : symmetry_numbers) {

    if ((0 > n_sym) || (n_sym >= n_symmetries_)) {
      Log.err("Error building subgroup of PermutationGroup: "
              "invalid symmetry index");
    }
    subgroup_permutations.push_back(permutations_[n_sym]);
  }

  return PermutationGroup(subgroup_permutations);
}

bool PermutationGroup::operator==(PermutationGroup const &rhs) const {
  return (permutations_ == rhs.permutations_);
}

bool PermutationGroup::operator!=(PermutationGroup const &rhs) const {
  return !operator==(rhs);
}

PermutationGroup allowed_subgroup(PermutationGroup const &group,
                                  Representation const &irrep) {

  auto const &allowed_symmetries = irrep.allowed_symmetries();

  if (allowed_symmetries.size() > 0) {
    return group.subgroup(allowed_symmetries);
  } else {
    return group;
  }
}

} // namespace hydra
