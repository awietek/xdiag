#pragma once

#include <hydra/symmetries/permutation_group.h>

namespace hydra::symmetries {

bool is_sublattice_stable(int n_sublat, PermutationGroup const &group);

std::vector<PermutationGroup>
get_sublattice_symmetries(int n_sublat, PermutationGroup const &group);

} // namespace hydra::symmetries
