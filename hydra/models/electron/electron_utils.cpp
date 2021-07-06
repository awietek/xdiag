#include "electron_utils.h"

#include <lila/utils/logger.h>

#include <hydra/common.h>

namespace hydra::electron {

void check_nup_ndn(int n_sites, int nup, int ndn) {
  if ((nup < 0) || (nup > n_sites))
    lila::Log.err("Error creating Electron: "
                  "invalid value of nup");
  if ((ndn < 0) || (ndn > n_sites))
    lila::Log.err("Error creating Electron: "
                  "invalid value of ndn");
}

} // namespace hydra::electron
