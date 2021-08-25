#pragma once

#include <string>
#include <vector>

namespace hydra::utils {

std::vector<std::vector<int>> read_permutations(std::string filename);

bool is_valid_permutation(int n_sites, const int *permutation);

template <class bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation);

template <class bit_t, class GroupAction>
std::vector<int> stabilizer_symmetries(bit_t bits, GroupAction &&group) {
  std::vector<int> stable_syms;
  for (int sym = 0; sym < group.n_symmetries(); ++sym)
    if (group.apply(sym, bits) == bits)
      stable_syms.push_back(sym);
  return stable_syms;
}

} // namespace hydra::utils
