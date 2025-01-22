#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag::testcases::tj {

OpSum tJchain(int nsites, double t, double J);
std::tuple<OpSum, arma::vec> tJchain_fullspectrum_alps(int L);
std::tuple<OpSum, arma::vec> tj_square2x2_fullspectrum_alps();
OpSum tj_alltoall(int nsites);
OpSum tj_alltoall_complex(int nsites);
std::tuple<OpSum, arma::vec> randomAlltoAll3();
std::tuple<OpSum, arma::vec> randomAlltoAll4();

} // namespace xdiag::testcases::tj
