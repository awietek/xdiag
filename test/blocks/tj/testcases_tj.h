#pragma once

#include <hydra/all.h>
#include <lila/all.h>

namespace hydra::testcases::tj {

std::tuple<BondList, Couplings> tJchain(int n_sites, double t, double J);

std::tuple<BondList, Couplings, lila::Vector<double>>
tJchain_fullspectrum_alps(int L);

std::tuple<BondList, Couplings, lila::Vector<double>>
tj_square2x2_fullspectrum_alps();

std::tuple<BondList, Couplings> tj_alltoall(int n_sites);
std::tuple<BondList, Couplings> tj_alltoall_complex(int n_sites);

std::tuple<BondList, Couplings, lila::Vector<double>> randomAlltoAll3();
std::tuple<BondList, Couplings, lila::Vector<double>> randomAlltoAll4();

} // namespace hydra::testcases::tj
