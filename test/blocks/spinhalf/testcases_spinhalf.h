#pragma once

#ifdef HYDRA_ENABLE_MPI
#include <mpi.h>
#endif

#include <hydra/all.h>
#include <lila/all.h>

namespace hydra::testcases::spinhalf {

std::tuple<BondList, Couplings> HBchain(int n_sites, double J1, double J2 = 0);

std::tuple<BondList, Couplings, lila::Vector<double>>
HBchain_fullspectrum_nup(int n_sites, int nup);

std::tuple<BondList, Couplings> HB_alltoall(int n_sites);

std::tuple<BondList, Couplings, double> triangular_12_complex(int nup,
                                                              double eta);

} // namespace hydra::testcases::spinhalf
