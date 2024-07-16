#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag::testcases::spinhalf {

OpSum HBchain(int64_t n_sites, double J1, double J2 = 0);
std::tuple<OpSum, arma::vec> HBchain_fullspectrum_nup(int64_t n_sites,
                                                      int64_t nup);
OpSum HB_alltoall(int64_t n_sites);
std::tuple<OpSum, double> triangular_12_complex(int64_t nup, double eta);

} // namespace xdiag::testcases::spinhalf
