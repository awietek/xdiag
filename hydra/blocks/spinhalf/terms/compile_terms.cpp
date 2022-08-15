#include "compile_terms.h"
#include <hydra/utils/logger.h>

namespace hydra::terms::spinhalf {

std::pair<BondList, Couplings> compile_terms(BondList const &bonds,
                                             Couplings const &couplings) {

  BondList bonds_compiled;
  Couplings couplings_compiled = couplings;

  for (auto bond : bonds) {

    auto type = bond.type();
    auto cpl = bond.coupling();
    auto sites = bond.sites();

    // Ignore bonds with zero coupling
    if (lila::close(couplings[cpl], 0.)) {
      continue;
    }

    // expand Heisenberg bonds into Ising and exchange
    if ((type == "HB") || (type == "HEISENBERG")) {
      bonds_compiled << Bond("ISING", cpl, sites);
      bonds_compiled << Bond("EXCHANGE", cpl, sites);
    } else if ((type == "ISING") || (type == "EXCHANGE") ||
               (type == "SCALARCHIRALITY") || (type == "S+") ||
               (type == "S-") || (type == "SZ")) {
      bonds_compiled << bond;
    } else {
      Log.err("Error in spinhalf::compile_terms: Unknown bond type {}",
              bond.type());
    }
  }
  return {bonds_compiled, couplings_compiled};
}

} // namespace hydra::terms::spinhalf
