#include "symmetric_operator.h"

namespace hydra {

std::pair<BondList, Couplings>
SymmetricOperator(BondList const &bonds, Couplings const &cpls,
                  PermutationGroup const &group) {
  BondList bonds_sym;
  Couplings cpls_sym;
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
      bonds_sym << Bond(type, cpl, sites_sym);
    }
    cpls_sym[cpl] = cpls[cpl] / (complex)N_group;
  }

  return {bonds_sym, cpls_sym};
}

} // namespace hydra
