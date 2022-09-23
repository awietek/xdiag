#include "symmetric_operator.h"

namespace hydra {

BondList symmetric_operator(BondList const &bonds,
                            PermutationGroup const &group) {
  BondList bonds_sym;
  int N_group = group.size();

  for (auto bond : bonds) {

    auto type = bond.type();
    auto cpl = bond.coupling();

    // Create all symmetrized bonds
    for (auto const &perm : group) {
      std::vector<int> sites_sym(bond.size(), 0);
      for (int site_idx = 0; site_idx < bond.size(); ++site_idx) {
        sites_sym[site_idx] = perm[bond[site_idx]];
      }
      bonds_sym << Bond(type, cpl / (complex)N_group, sites_sym);
    }
  }

  return bonds_sym;
}

} // namespace hydra
