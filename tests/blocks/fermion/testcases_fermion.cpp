// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <random>

#include "testcases_fermion.hpp"

namespace xdiag::testcases::fermion {

OpSum freefermion_alltoall(int64_t nsites) {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0., 1.);
  OpSum ops;
  for (int64_t s1 = 0; s1 < nsites; ++s1)
    for (int64_t s2 = s1 + 1; s2 < nsites; ++s2) {
      double c = distribution(generator);
      ops += c * Op("Hop", {s1, s2});
    }
  return ops;
}

OpSum freefermion_alltoall_complex(int64_t nsites) {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0., 1.);
  OpSum ops;
  for (int64_t s1 = 0; s1 < nsites; ++s1)
    for (int64_t s2 = s1 + 1; s2 < nsites; ++s2) {
      complex c = complex(distribution(generator), distribution(generator));
      ops += std::real(c) * Op("Hop", {s1, s2});
      ops += std::complex(0.0, std::imag(c)) * Op("HopAsym", {s1, s2});
    }
  return ops;
}

} // namespace xdiag::testcases::fermion
