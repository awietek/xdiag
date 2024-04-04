#pragma once

#include <xdiag/operators/bondlist.hpp>

namespace xdiag::testcases::tj {

BondList tJchain(int n_sites, double t, double J);
std::tuple<BondList, arma::vec> tJchain_fullspectrum_alps(int L);
std::tuple<BondList, arma::vec> tj_square2x2_fullspectrum_alps();
BondList tj_alltoall(int n_sites);
BondList tj_alltoall_complex(int n_sites);
std::tuple<BondList, arma::vec> randomAlltoAll3();
std::tuple<BondList, arma::vec> randomAlltoAll4();

} // namespace xdiag::testcases::tj
