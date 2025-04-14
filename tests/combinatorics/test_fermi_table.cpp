// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <iostream>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <xdiag/combinatorics/fermi_table.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/symmetries/operations/fermi_sign.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

template <typename bit_t>
void test_fermi_bool_table(PermutationGroup const &group) {
  using combinatorics::Combinations;
  using combinatorics::Subsets;
  using namespace symmetries;

  int nsites = group.nsites();
  int n_symmetries = group.size();

  for (int npar = 0; npar <= nsites; ++npar) {

    auto fermi_tbl =
        combinatorics::FermiTableCombinations<bit_t>(nsites, npar, group);
    for (int sym = 0; sym < n_symmetries; ++sym) {
      for (bit_t state : Combinations<bit_t>(nsites, npar)) {
        REQUIRE(fermi_tbl.sign(sym, state) ==
                fermi_bool_of_permutation(state, group[sym]));
      }
    }
  }

  auto fermi_tbl = combinatorics::FermiTableSubsets<bit_t>(nsites, group);
  for (int sym = 0; sym < n_symmetries; ++sym) {
    for (bit_t state : Subsets<bit_t>(nsites)) {
      REQUIRE(fermi_tbl.sign(sym, state) ==
              fermi_bool_of_permutation(state, group[sym]));
    }
  }
}
TEST_CASE("fermi_table", "[symmetries]") {
  xdiag::Log("Test fermi_table");
  int max_N = 6;

  for (int nsites = 1; nsites <= max_N; ++nsites) {
    Log("chain N={}", nsites);
    auto irreps = xdiag::testcases::electron::get_cyclic_group_irreps(nsites);
    test_fermi_bool_table<uint16_t>(irreps[0].group());
    test_fermi_bool_table<uint32_t>(irreps[0].group());
    test_fermi_bool_table<uint64_t>(irreps[0].group());
  }

  Log("triangular 3x3");
  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";
  auto fl = FileToml(lfile);
  auto group = fl["Symmetries"].as<PermutationGroup>();
  test_fermi_bool_table<uint16_t>(group);
  test_fermi_bool_table<uint32_t>(group);
  test_fermi_bool_table<uint64_t>(group);

  Log("done");
}
