// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

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

std::pair<int, int> target_nup_ndn(std::string op_str, int nup, int ndn);
std::pair<int, int> target_nup_ndn(std::string op_str1, std::string op_str2,
                                   int nup, int ndn);

bool valid_nup_ndn(int nup, int ndn, int nsites);
bool valid_nup_ndn(std::string op_str, int nup, int ndn, int nsites);
bool valid_nup_ndn(std::string op_str1, std::string op_str2, int nup, int ndn,
                   int nsites);

} // namespace xdiag::testcases::tj
