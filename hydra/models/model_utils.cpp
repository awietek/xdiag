#include "model_utils.h"

namespace hydra::utils {

bool coupling_is_zero(Bond const &bond, Couplings const &couplings) {
  std::string coupling = bond.coupling();
  return (!couplings.defined(coupling)) ||
         lila::close(couplings[coupling], (complex)0.);
}

void check_nup_ndn_electron(int n_sites, int nup, int ndn) {
  if ((nup < 0) || (nup > n_sites))
    lila::Log.err("Error creating Electron: "
                  "invalid value of nup");
  if ((ndn < 0) || (ndn > n_sites))
    lila::Log.err("Error creating Electron: "
                  "invalid value of ndn");
}

void check_nup_ndn_tj(int n_sites, int nup, int ndn) {
  if ((nup < 0) || (nup > n_sites))
    lila::Log.err("Error creating tJ: "
                  "invalid value of nup");
  if ((ndn < 0) || (ndn > n_sites))
    lila::Log.err("Error creating tJ: "
                  "invalid value of ndn");
  if (nup + ndn > n_sites)
    lila::Log.err("Error creating tJ: "
                  "invalid value of nup+ndn");
}

} // namespace hydra::utils
