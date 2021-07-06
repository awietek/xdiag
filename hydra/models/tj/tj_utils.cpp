#include "tj_utils.h"

#include <lila/utils/logger.h>

#include <hydra/common.h>

namespace hydra::tjdetail {

void check_nup_ndn(int n_sites, int nup, int ndn) {
  if ((nup < 0) || (nup > n_sites))
    lila::Log.err("Error creating tJ: "
                 "invalid value of nup");
  if ((ndn < 0) || (ndn > n_sites))
    lila::Log.err("Error creating tJ: "
                 "invalid value of ndn");
  if (nup+ndn > n_sites)
    lila::Log.err("Error creating tJ: "
                 "nup+ndn > n_sites");
}



} // namespace hydra::tjdetail
