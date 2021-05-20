#include "testcases_tj.h"

namespace hydra::testcases::tj {

std::tuple<BondList, Couplings> tJchain(int n_sites, double t, double J) {

  BondList bondlist;
  Couplings couplings;
  couplings["T"] = t;
  couplings["J"] = J;
  for (int s = 0; s < n_sites; ++s) {
    bondlist << Bond("HOP", "T", {s, (s + 1) % n_sites});
    bondlist << Bond("HB", "J", {s, (s + 1) % n_sites});
  }
  return std::make_tuple(bondlist, couplings);
}

} // namespace hydra::testcases::tj
