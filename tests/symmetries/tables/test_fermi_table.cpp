// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <iostream>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/catch.hpp>

#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/config.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/symmetries/fermi_sign.hpp>
#include <xdiag/symmetries/tables/fermi_table.hpp>
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
        symmetries::FermiTable(Combinations<bit_t>(nsites, npar), group);
    for (int sym = 0; sym < n_symmetries; ++sym) {
      for (bit_t state : Combinations<bit_t>(nsites, npar)) {
        REQUIRE(fermi_tbl.sign(sym, state) ==
                fermi_bool_of_permutation(state, group[sym]));
      }
    }
  }

  auto fermi_tbl = symmetries::FermiTable(Subsets<bit_t>(nsites), group);
  for (int sym = 0; sym < n_symmetries; ++sym) {
    for (bit_t state : Subsets<bit_t>(nsites)) {
      REQUIRE(fermi_tbl.sign(sym, state) ==
              fermi_bool_of_permutation(state, group[sym]));
    }
  }
}
TEST_CASE("fermi_table", "[symmetries]") {
  Log("Test fermi_table");
  int max_N = 6;

  for (int nsites = 1; nsites <= max_N; ++nsites) {
    Log("  chain N={}", nsites);
    auto irreps = xdiag::testcases::electron::get_cyclic_group_irreps(nsites);
    test_fermi_bool_table<uint32_t>(irreps[0].group());
    test_fermi_bool_table<uint64_t>(irreps[0].group());
  }

  Log("  triangular 3x3");
  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";
  auto fl = FileToml(lfile);
  auto group = fl["Symmetries"].as<PermutationGroup>();
  test_fermi_bool_table<uint32_t>(group);
  test_fermi_bool_table<uint64_t>(group);

  Log("  long chain");
  {
    using namespace bits;
    using namespace combinatorics;
    using namespace symmetries;
    int nsites = 128;
    auto irreps = xdiag::testcases::electron::get_cyclic_group_irreps(nsites);
    auto group = irreps[0].group();
    int n_symmetries = group.size();
    for (int npar = 0; npar <= 2; ++npar) {

      auto fermi_tbl =
          FermiTable(Combinations<BitsetDynamic>(nsites, npar), group);
      for (int sym = 0; sym < n_symmetries; ++sym) {
        for (auto state : Combinations<BitsetDynamic>(nsites, npar)) {
          REQUIRE(fermi_tbl.sign(sym, state) ==
                  fermi_bool_of_permutation(state, group[sym]));
        }
      }
    }
  }
}
