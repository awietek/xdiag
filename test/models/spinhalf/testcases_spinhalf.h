#pragma once

#include <hydra/all.h>
#include <lila/all.h>

namespace hydra::testcases::spinhalf {

std::tuple<BondList, Couplings> HBchain(int n_sites, double J);

std::tuple<BondList, Couplings, lila::Vector<double>>
HBchain_fullspectrum_nup(int n_sites, int nup);

std::tuple<BondList, Couplings> HB_alltoall(int n_sites);

  
} // namespace hydra::testcases::spinhalf
