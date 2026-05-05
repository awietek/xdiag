// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <tuple>

#include <xdiag/armadillo.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::testcases::spinhalf {

// Heisenberg chain with nearest-neighbor J1 and next-nearest-neighbor J2
OpSum HBchain(int64_t nsites, double J1, double J2 = 0);

// Heisenberg chain J1=1, J2=0 with exact full spectrum for given nup
std::tuple<OpSum, arma::vec> HBchain_fullspectrum_nup(int64_t nsites,
                                                      int64_t nup);

// Heisenberg all-to-all with random couplings
OpSum HB_alltoall(int64_t nsites);

// Triangular N=12 lattice with complex exchange; returns ops and exact e0
std::tuple<OpSum, double> triangular_12_complex(int64_t nup, double eta);

} // namespace xdiag::testcases::spinhalf
