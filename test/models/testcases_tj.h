#pragma once

#include <hydra/all.h>

namespace hydra::testcases::tj {

std::tuple<BondList, Couplings> tJchain(int n_sites, double t, double J);

} // namespace hydra::testcases::tj
