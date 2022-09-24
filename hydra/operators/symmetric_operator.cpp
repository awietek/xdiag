#include "symmetric_operator.h"

namespace hydra {

BondList symmetric_operator(BondList const &bonds,
                            PermutationGroup const &group) {
  BondList bonds_sym;
  int N_group = group.size();

  for (auto bond : bonds) {

    complex coupling = 0.;

    if (bond.coupling_named()) {
      std::string name = bond.coupling_name();
      if (bonds.coupling_defined(name)) {
        coupling = bonds.coupling(name);
      } else {
        Log.err("Error in symmetric_operator: coupling with name {} undefined "
                "in BondList",
                name);
      }
    } else {
      coupling = bond.coupling();
    }

    if (bond.type_defined()) {
      std::string type = bond.type();

      // Create all symmetrized bonds
      for (auto const &perm : group) {
        std::vector<int> sites_sym(bond.size(), 0);
        for (int site_idx = 0; site_idx < bond.size(); ++site_idx) {
          sites_sym[site_idx] = perm[bond[site_idx]];
        }
        bonds_sym << Bond(type, coupling / (complex)N_group, sites_sym);
      }
    } else {
      arma::cx_mat mat = bond.matrix();

      // Create all symmetrized bonds
      for (auto const &perm : group) {
        std::vector<int> sites_sym(bond.size(), 0);
        for (int site_idx = 0; site_idx < bond.size(); ++site_idx) {
          sites_sym[site_idx] = perm[bond[site_idx]];
        }
        bonds_sym << Bond(mat, coupling / (complex)N_group, sites_sym);
      }
    }
  }

  return bonds_sym;
}

} // namespace hydra
