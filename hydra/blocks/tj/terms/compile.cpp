#include "compile.h"
#include <hydra/operators/compiler.h>
#include <hydra/utils/logger.h>
#include <hydra/utils/print_macro.h>

namespace hydra::tj {

BondList compile(BondList const &bonds, double precision) {
  using namespace operators;

  BondList bonds_explicit = compile_explicit(bonds, precision, "keep");

  BondList bonds_compiled;
  for (auto bond : bonds_explicit) {
    if (bond.type_defined()) {
      std::string type = bond.type();
      if (std::find(special_bond_types.begin(), special_bond_types.end(),
                    type) == special_bond_types.end()) {
        Log.err("Error compiling BondList: invalid or undefined type found: {}",
                type);
      } else {

        // Exchange and Ising terms
        if ((type == "HB") || (type == "HEISENBERG")) {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("ISING", bond.coupling(), bond.sites());
          bonds_compiled << Bond("EXCHANGE", bond.coupling(), bond.sites());
        } else if ((type == "TJHB") || (type == "TJHEISENBERG")) {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("TJISING", bond.coupling(), bond.sites());
          bonds_compiled << Bond("EXCHANGE", bond.coupling(), bond.sites());
        } else if (type == "ISING") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("ISING", bond.coupling(), bond.sites());
        } else if (type == "EXCHANGE") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("EXCHANGE", bond.coupling(), bond.sites());
        }

        // Hopping terms
        else if (type == "HOP") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("HOPUP", bond.coupling(), bond.sites());
          bonds_compiled << Bond("HOPDN", bond.coupling(), bond.sites());
        } else if (type == "HOPUP") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("HOPUP", bond.coupling(), bond.sites());
        } else if (type == "HOPDN") {
          check_bond_has_correct_number_of_sites(bond, 2);
          check_bond_has_disjoint_sites(bond);
          bonds_compiled << Bond("HOPDN", bond.coupling(), bond.sites());
        }

        // Number operators
        else if (type == "NUMBER") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("NUMBERUP", bond.coupling(), bond.sites());
          bonds_compiled << Bond("NUMBERDN", bond.coupling(), bond.sites());
        } else if (type == "NUMBERUP") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("NUMBERUP", bond.coupling(), bond.sites());
        } else if (type == "NUMBERDN") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("NUMBERDN", bond.coupling(), bond.sites());
        } else if (type == "SZ") {
          check_bond_has_correct_number_of_sites(bond, 1);
          bonds_compiled << Bond("NUMBERUP", 0.5 * bond.coupling(),
                                 bond.sites());
          bonds_compiled << Bond("NUMBERDN", -0.5 * bond.coupling(),
                                 bond.sites());
        }
      }
    }
  }
  return bonds_compiled;
}

} // namespace hydra::tj
