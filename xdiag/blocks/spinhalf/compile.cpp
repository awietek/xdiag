#include "compile.hpp"

#include <xdiag/operators/compiler.hpp>
#include <xdiag/operators/non_branching_bonds.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/print_macro.hpp>
#include <xdiag/utils/timing.hpp>

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
        XDIAG_THROW(std::string("Invalid or undefined type: \"") + type +
                    std::string("\""));
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
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return BondList();
}

} // namespace xdiag::spinhalf
