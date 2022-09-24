#include "compile.h"
#include <hydra/operators/compiler.h>
#include <hydra/utils/logger.h>

namespace hydra::spinhalf {

BondList compile(BondList const &bonds, double precision) {

  BondList bonds_explicit =
      operators::compile_explicit(bonds, precision, "keep");
  BondList bonds_special;
  BondList bonds_generic;
  for (auto bond : bonds_explicit) {
    if (bond.type_defined()) {
      std::string type = bond.type();
      if (std::find(special_bond_types.begin(), special_bond_types.end(),
                    type) == special_bond_types.end()) {
        Log.err("Error compiling BondList: invalid or undefined type found: {}",
                type);
      } else {

        if ((type == "HB") || (type == "HEISENBERG")) {
          bonds_special << Bond("ISING", bond.coupling(), bond.sites());
          bonds_special << Bond("EXCHANGE", bond.coupling(), bond.sites());
        } else {
          bonds_special << bond;
        }
      }
    } else {
      bonds_generic << bond;
    }
  }
  return bonds_special + bonds_generic;
}

} // namespace hydra::spinhalf
