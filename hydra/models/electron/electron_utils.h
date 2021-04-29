#pragma once

#include <vector>

namespace hydra::electron {
void check_nup_ndn(int n_sites, int nup, int ndn);

template <class bit_t, class SymmetryGroup>
std::vector<int> stabilizer_symmetries(bit_t bits, SymmetryGroup &&group) {
  std::vector<int> stable_syms;
  for (int sym = 0; sym < group.size(); ++sym)
    if (group.apply(sym, bits) == bits)
      stable_syms.push_back(sym);
  return stable_syms;
}

} // namespace hydra::electron
