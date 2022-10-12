#include "symmetrized_operator.h"

namespace hydra {

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
                              Representation const &irrep) {
  BondList bonds_sym;
  int N_group = group.size();

  for (auto bond : bonds) {

    complex coupling = 0.;

    if (bond.coupling_named()) {
      std::string name = bond.coupling_name();
      if (bonds.coupling_defined(name)) {
        coupling = bonds.coupling(name);
      } else {
        Log.err(
            "Error in symmetrized_operator: coupling with name {} undefined "
            "in BondList",
            name);
      }
    } else {
      coupling = bond.coupling();
    }

    if (bond.type_defined()) {
      std::string type = bond.type();

      // Create all symmetrized bonds
      for (int i = 0; i < N_group; ++i) {
        Permutation perm = group[i];
        complex bloch = irrep.character(i);

        std::vector<int> sites_sym(bond.size(), 0);
        for (int site_idx = 0; site_idx < bond.size(); ++site_idx) {
          sites_sym[site_idx] = perm[bond[site_idx]];
        }
        complex cpl_sym = bloch * coupling / (complex)N_group;
        bonds_sym << Bond(type, cpl_sym, sites_sym);
      }
    } else {
      arma::cx_mat mat = bond.matrix();

      // Create all symmetrized bonds
      for (int i = 0; i < N_group; ++i) {
        Permutation perm = group[i];
        complex bloch = irrep.character(i);

        std::vector<int> sites_sym(bond.size(), 0);
        for (int site_idx = 0; site_idx < bond.size(); ++site_idx) {
          sites_sym[site_idx] = perm[bond[site_idx]];
        }

        complex cpl_sym = bloch * coupling / (complex)N_group;
        bonds_sym << Bond(mat, cpl_sym, sites_sym);
      }
    }
  }

  return bonds_sym;
}

} // namespace hydra
