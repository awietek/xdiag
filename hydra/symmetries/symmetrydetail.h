#pragma once

#include <vector>

namespace hydra {
namespace detail {

bool is_valid_permutation(std::vector<int> const& permutation);

template <class bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation) {
  bit_t tstate = 0;
  for (int i = 0; i < n_sites; ++i) {
    int val = siteval(state, i);
    set_siteval(&tstate, permutation[i], val);
  }
  return tstate;
}

double fermi_sign(bit_t state, int n_sites, const int *permutation);

} // namespace detail
} // namespace hydra
