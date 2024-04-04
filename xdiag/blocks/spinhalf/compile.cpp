#include "compile.h"
#include <xdiag/operators/compiler.h>
#include <xdiag/operators/non_branching_bonds.h>
#include <xdiag/utils/logger.h>
#include <xdiag/utils/print_macro.h>
#include <xdiag/utils/timing.h>

namespace xdiag::spinhalf {

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
        XDiagThrow(std::runtime_error,
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
  XDiagRethrow("Unable to compile BondList");
  return BondList();
}

} // namespace xdiag::spinhalf
