#include "symmetrized_operator.hpp"

#include <xdiag/operators/compiler.hpp>

namespace xdiag {

BondList symmetrized_operator(Bond const &bond, PermutationGroup const &group) {
  auto bonds = BondList({bond});
  return symmetrized_operator(bonds, group);
}

BondList symmetrized_operator(Bond const &bond, PermutationGroup const &group,
                              Representation const &irrep) {
  auto bonds = BondList({bond});
  return symmetrized_operator(bonds, group, irrep);
}

BondList symmetrized_operator(BondList const &bonds,
                              PermutationGroup const &group) {
  return symmetrized_operator(bonds, group,
                              trivial_representation(group.size()));
}

BondList symmetrized_operator(BondList const &bonds,
                              PermutationGroup const &group,
                              Representation const &irrep) try {
  BondList bonds_sym;
  int64_t N_group = group.size();
  BondList bonds_explicit = make_explicit(bonds);
  for (auto bond : bonds_explicit) {

    std::string type = bond.type();
    Coupling coupling = bond.coupling();

    // Create all symmetrized bonds
    for (int64_t i = 0; i < N_group; ++i) {
      Permutation perm = group[i];
      complex bloch = irrep.character(i);

      std::vector<int64_t> sites_sym(bond.size(), 0);
      for (int64_t site_idx = 0; site_idx < bond.size(); ++site_idx) {
        sites_sym[site_idx] = perm[bond[site_idx]];
      }
      Coupling cpl_sym = bloch * coupling / (complex)N_group;
      bonds_sym += Bond(type, cpl_sym, sites_sym);
    }
  }
  return bonds_sym;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
