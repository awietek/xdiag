#include "compile.h"
#include <hydra/operators/compiler.h>
#include <hydra/utils/logger.h>

namespace hydra::electron {

BondList compile(BondList const &bonds, double precision) {

  BondList bonds_explicit =
      operators::compile_explicit(bonds, precision, "keep");
  
  BondList bonds_special;
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
        } else if (type == "HOP") {
          bonds_special << Bond("HOPUP", bond.coupling(), bond.sites());
          bonds_special << Bond("HOPDN", bond.coupling(), bond.sites());
        } else {
          bonds_special << bond;
        }
      }
    }
  }

  // Set Hubbbbard U term again
  if (bonds.coupling_defined("U")) {
    complex coupling = bonds.coupling("U");
    if (std::abs(coupling) > precision) {
      bonds_special["U"] = coupling;
    }
  }

  return bonds_special;
}

} // namespace hydra::electron
