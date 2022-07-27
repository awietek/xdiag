#pragma once

#include <hydra/symmetries/permutation_group.h>
#include <lila/external/gsl/span>

#define HYDRA_SUBLATTICE_UNSTABLE -1

namespace hydra::symmetries {

bool is_sublattice_permutation(int n_sublat, int sublat,
                               gsl::span<int const> permutation);

int which_sublattice_permutation(int n_sublat,
                                 gsl::span<int const> permutation);

bool is_sublattice_stable(int n_sublat, PermutationGroup const &group);

std::vector<int> sublattice_permutations(int n_sublat, int sublat,
                                         PermutationGroup const &group);

} // namespace hydra::symmetries