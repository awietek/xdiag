// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include <xdiag/blocks/tj.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>

using namespace xdiag;

void check_dimensions_sum_up_tj_symmetric(int64_t nsites,
                                          std::vector<Representation> irreps) {
  using combinatorics::binomial;

  Log.out("tj_symmetric: dimension sum test. N: {}", nsites);

  int64_t sum_of_dims = 0;

  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn <= nsites - nup; ++ndn) {
      int64_t sum_of_dims_updn = 0;
      for (auto irrep : irreps) {

        auto block = tJ(nsites, nup, ndn, irrep);
        sum_of_dims += block.size();
        sum_of_dims_updn += block.size();
      }
      REQUIRE(sum_of_dims_updn ==
              binomial(nsites, nup) * binomial(nsites - nup, ndn));
    }
  }
  int64_t p = 1;
  for (int64_t i = 0; i < nsites; ++i) {
    p *= 3;
  }
  REQUIRE(sum_of_dims == p);
}

TEST_CASE("tj_symmetric", "[tj]") {

  // Test linear chains
  Log("tj_symmetric: chain test");
  for (int64_t nsites = 1; nsites < 7; ++nsites) {

    // test cyclic group
    std::vector<Permutation> permutation_array;
    for (int64_t sym = 0; sym < nsites; ++sym) {

      std::vector<int64_t> pv;
      for (int64_t site = 0; site < nsites; ++site) {
        int64_t newsite = (site + sym) % nsites;
        pv.push_back(newsite);
      }
      permutation_array.push_back(Permutation(pv));
    }
    auto group = PermutationGroup(permutation_array);

    // Create irreps
    std::vector<Representation> irreps;
    for (int64_t k = 0; k < nsites; ++k) {
      std::vector<complex> chis;
      for (int64_t l = 0; l < nsites; ++l)
        chis.push_back({std::cos(2 * M_PI * l * k / nsites),
                        std::sin(2 * M_PI * l * k / nsites)});
      irreps.push_back(Representation(group, chis));
    }

    check_dimensions_sum_up_tj_symmetric(nsites, irreps);
  }

  // test a 3x3 triangular lattice
  int64_t nsites = 9;
  Log("tj_symmetric: triangular 3x3 test");

  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";
  auto fl = FileToml(lfile);
  std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};

  std::vector<Representation> irreps;
  for (auto [name, mult] : rep_name_mult) {
    for (int64_t i = 0; i < mult; ++i) {
      irreps.push_back(read_representation(fl, name));
    }
  }

  check_dimensions_sum_up_tj_symmetric(nsites, irreps);
}
