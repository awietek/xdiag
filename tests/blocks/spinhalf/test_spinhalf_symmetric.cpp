// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>

using namespace xdiag;
using namespace xdiag::combinatorics;

void test_spinchain_blocks(int64_t nsites) {

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

  // Create irrep with K momentum
  int64_t sum_of_dims = 0;
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    int64_t sum_of_dims_up = 0;

    for (int64_t k = 0; k < nsites; ++k) {

      // Create characters
      std::vector<complex> chis;
      for (int64_t l = 0; l < nsites; ++l)
        chis.push_back({std::cos(2 * M_PI * l * k / nsites),
                        std::sin(2 * M_PI * l * k / nsites)});
      auto irrep = Representation(group, chis);

      auto block = Spinhalf(nsites, nup, irrep);

      sum_of_dims += block.size();
      sum_of_dims_up += block.size();
    }
    REQUIRE(sum_of_dims_up == binomial(nsites, nup));
  }
  REQUIRE(sum_of_dims == pow(2, nsites));
}

TEST_CASE("spinhalf_symmetric", "[spinhalf]") {

  // Test the tJ chain
  for (int64_t nsites = 1; nsites < 8; ++nsites) {
    Log.out("spinhalf_symmetric: block test: Spinhalf Chain {}", nsites);
    test_spinchain_blocks(nsites);
  }

  // test a 3x3 triangular lattice
  Log.out("spinhalf_symmetric: block test: Triangular 3x3");
  int64_t nsites = 9;

  std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
      {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
      {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
      {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
      {"Y.D1.A", 6},      {"Y.D1.B", 6}};

  std::string lfile = XDIAG_DIRECTORY
      "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
  auto fl = FileToml(lfile);
  int64_t sum_dim = 0;
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    int64_t sum_dim_up = 0;
    for (auto [name, mult] : rep_name_mult) {
      auto irrep = read_representation(fl, name);
      auto block = Spinhalf(nsites, nup, irrep);
      int64_t dim = block.size() * mult;
      sum_dim_up += dim;
      sum_dim += dim;
    }
    REQUIRE(sum_dim_up == binomial(nsites, nup));
  }
  REQUIRE(sum_dim == (int64_t)pow(2, nsites));
}
