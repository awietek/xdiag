// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sublattice_stability.hpp"

#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::symmetries {

bool is_sublattice_permutation(int n_sublat, int sublat,
                               Permutation const &permutation) {
  int nsites = permutation.size();
  int nsites_sublat = nsites / n_sublat;

  // permutation is sublattice stable, if the sublattice
  // (i.e. indices between begin and end) are mapped to
  // most significant indices (nsites-nsites_sublat) to nsites
  int begin = sublat * nsites_sublat;
  int end = (sublat + 1) * nsites_sublat;
  for (int j = begin; j < end; ++j) {
    if (permutation[j] < nsites - nsites_sublat) {
      return false;
    }
  }
  return true;
}

int which_sublattice_permutation(int n_sublat, Permutation const &permutation) {
  int sublat = XDIAG_SUBLATTICE_UNSTABLE;
  for (int s = 0; s < n_sublat; ++s) {
    if (is_sublattice_permutation(n_sublat, s, permutation)) {
      sublat = s;
      break;
    }
  }
  return sublat;
}

bool is_sublattice_stable(int n_sublat, PermutationGroup const &group) {

  for (auto const &perm : group) {
    int sublat = which_sublattice_permutation(n_sublat, perm);
    if (sublat == XDIAG_SUBLATTICE_UNSTABLE) {
      return false;
    }
  }
  return true;
}

std::vector<int> sublattice_permutations(int n_sublat, int sublat,
                                         PermutationGroup const &group) {
  std::vector<int> syms;

  int sym = 0;
  for (auto const &perm : group) {
    if (which_sublattice_permutation(n_sublat, perm) == sublat) {
      syms.push_back(sym);
    }
    ++sym;
  }
  return syms;
}

} // namespace xdiag::symmetries
