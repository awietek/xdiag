#pragma once

#include <vector>
#include <string>

namespace hydra {
std::vector<std::vector<int>> read_permutations(std::string filename);

namespace detail {

bool is_valid_permutation(int n_sites, const int *permutation);

template <class bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation);
  
template <class bit_t>
double fermi_sign(bit_t state, int n_sites, const int *permutation);

} // namespace detail
} // namespace hydra
