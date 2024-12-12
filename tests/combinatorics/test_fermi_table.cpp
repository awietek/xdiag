#include "../catch.hpp"

#include <iostream>

#include <tests/blocks/electron/testcases_electron.hpp>
#include <xdiag/combinatorics/fermi_table.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/symmetries/operations/fermi_sign.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>

using namespace xdiag;

template <typename bit_t>
void test_fermi_bool_table(PermutationGroup const &group) {
  using combinatorics::Combinations;
  using combinatorics::Subsets;
  using namespace symmetries;

  int n_sites = group.n_sites();
  int n_symmetries = group.n_symmetries();

  for (int npar = 0; npar <= n_sites; ++npar) {

    auto fermi_tbl =
        combinatorics::FermiTableCombinations<bit_t>(n_sites, npar, group);
    for (int sym = 0; sym < n_symmetries; ++sym) {
      for (bit_t state : Combinations<bit_t>(n_sites, npar)) {
        REQUIRE(fermi_tbl.sign(sym, state) ==
                fermi_bool_of_permutation(state, group[sym]));
      }
    }
  }

  auto fermi_tbl = combinatorics::FermiTableSubsets<bit_t>(n_sites, group);
  for (int sym = 0; sym < n_symmetries; ++sym) {
    for (bit_t state : Subsets<bit_t>(n_sites)) {
      REQUIRE(fermi_tbl.sign(sym, state) ==
              fermi_bool_of_permutation(state, group[sym]));
    }
  }
}
TEST_CASE("fermi_table", "[symmetries]") {
  xdiag::Log("Test fermi_table");
  int max_N = 6;

  for (int n_sites = 0; n_sites <= max_N; ++n_sites) {
    Log("chain N={}", n_sites);
    auto [group, irreps] =
        xdiag::testcases::electron::get_cyclic_group_irreps(n_sites);
    (void)irreps;
    test_fermi_bool_table<uint16_t>(group);
    test_fermi_bool_table<uint32_t>(group);
    test_fermi_bool_table<uint64_t>(group);
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
