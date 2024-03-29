#include "compile.h"
#include <hydra/operators/compiler.h>
#include <hydra/operators/non_branching_bonds.h>
#include <hydra/utils/logger.h>
#include <hydra/utils/print_macro.h>
#include <hydra/utils/timing.h>

namespace hydra::spinhalf {

BondList compile(BondList const &bonds, double precision) try {

  BondList bonds_explicit =
      operators::compile_explicit(bonds, precision, "keep");
  BondList bonds_special;
  BondList bonds_generic;

  for (auto bond : bonds_explicit) {
    if (bond.type_defined()) {
      std::string type = bond.type();
      if (std::find(special_bond_types.begin(), special_bond_types.end(),
                    type) == special_bond_types.end()) {
        HydraThrow(std::runtime_error,
                   std::string("Invalid or undefined type found ") + type);
      } else {

        if ((type == "HB") || (type == "HEISENBERG")) {
          bonds_special << Bond("ISING", bond.coupling(), bond.sites());
          bonds_special << Bond("EXCHANGE", bond.coupling(), bond.sites());
        } else {
          bonds_special << bond;
        }
      }
    } else {
      BondList bonds_nb = operators::non_branching_bonds(bond, precision);
      for (auto b : bonds_nb) {
        bonds_generic << b;
      }
    }
  }
  return bonds_special + bonds_generic;
} catch (...) {
  HydraRethrow("Unable to compile BondList");
  return BondList();
}

} // namespace hydra::spinhalf
