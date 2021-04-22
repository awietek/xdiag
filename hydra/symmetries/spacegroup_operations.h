#pragma once

#include <string>
#include <vector>

namespace hydra {
std::vector<std::vector<int>> read_permutations(std::string filename);

namespace symmetries {

bool is_valid_permutation(int n_sites, const int *permutation);

template <class bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation);

} // namespace symmetries
} // namespace hydra
