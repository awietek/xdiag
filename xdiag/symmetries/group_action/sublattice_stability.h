#pragma once

#include <xdiag/extern/gsl/span>
#include <xdiag/symmetries/permutation_group.h>

#define XDIAG_SUBLATTICE_UNSTABLE -1

namespace xdiag::symmetries {

bool is_sublattice_permutation(int n_sublat, int sublat,
                               Permutation const &permutation);

int which_sublattice_permutation(int n_sublat, Permutation const &permutation);

bool is_sublattice_stable(int n_sublat, PermutationGroup const &group);

std::vector<int> sublattice_permutations(int n_sublat, int sublat,
                                         PermutationGroup const &group);

} // namespace xdiag::symmetries
