// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>
#include <xdiag/symmetries/group_action/group_action_sublattice.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

template <class Action1, class Action2>
void compare_actions(Action1 &&action1, Action2 &&action2) {
  int nsites = action1.nsites();
  int n_symmetries = action1.n_symmetries();

  REQUIRE(action2.nsites() == nsites);
  REQUIRE(action2.n_symmetries() == n_symmetries);

  for (auto bits : combinatorics::Subsets<uint64_t>(nsites)) {

    // Check translations
    for (int sym = 0; sym < n_symmetries; ++sym) {
      auto t1 = action1.apply(sym, bits);
      auto t2 = action2.apply(sym, bits);
      // Log("1sl bits: {} t1: {} t2: {}", BSTR(bits), BSTR(t1), BSTR(t2));
      REQUIRE(t1 == t2);
    }

    // Check representative search
    auto r1 = action1.representative(bits);
    auto r2 = action2.representative(bits);
    // Log("1sl bits: {} r1: {} r2: {}", BSTR(bits), BSTR(r1), BSTR(r2));
    REQUIRE(r1 == r2);

    {
      auto [r1, sym1] = action1.representative_sym(bits);
      auto [r2, sym2] = action2.representative_sym(bits);
      REQUIRE(r1 == r2);
      REQUIRE(action1.apply(sym1, bits) == r1);
      REQUIRE(action2.apply(sym2, bits) == r2);
    }

    {
      auto [r1, syms1] = action1.representative_syms(bits);
      auto [r2, syms2] = action2.representative_syms(bits);
      REQUIRE(r1 == r2);
      REQUIRE(syms1.size() == syms2.size());
      for (std::size_t i = 0; i < syms1.size(); ++i) {
        int sym1 = syms1[i];
        int sym2 = syms2[i];
        REQUIRE(action1.apply(sym1, bits) == r1);
        REQUIRE(action2.apply(sym2, bits) == r2);
      }
    }
  }
}

template <class bit_t> void test_group_action_sublattice() {

  {
    Log("GroupActionSublattice: 1 sublattice");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.2sl.toml";
    auto fl = FileToml(lfile);
    auto group = fl["Symmetries"].as<PermutationGroup>();
    auto action = GroupAction(group);
    auto action_sl = GroupActionSublattice<bit_t, 1>(group);
    compare_actions(action, action_sl);
  }

  // Two sublattice
  {
    Log("GroupActionSublattice: 2 sublattice");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.2sl.toml";
    auto fl = FileToml(lfile);
    auto group = fl["Symmetries"].as<PermutationGroup>();
    auto action = GroupAction(group);
    auto action_sl = GroupActionSublattice<bit_t, 2>(group);
    compare_actions(action, action_sl);
  }

  // Three sublattice
  {
    Log("GroupActionSublattice: 3 sublattice");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.9.heisenberg.3sl.toml";
    auto fl = FileToml(lfile);
    auto group = fl["Symmetries"].as<PermutationGroup>();
    auto action = GroupAction(group);
    auto action_sl = GroupActionSublattice<bit_t, 3>(group);
    compare_actions(action, action_sl);
  }

  // Three sublattice (triangular example)
  {
    Log("GroupActionSublattice: 3 sublattice (triangular)");
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
    auto fl = FileToml(lfile);
    auto group = fl["Symmetries"].as<PermutationGroup>();
    auto action = GroupAction(group);
    auto action_sl = GroupActionSublattice<bit_t, 3>(group);
    compare_actions(action, action_sl);
  }

  // Four sublattice
  {
    Log("GroupActionSublattice: 4 sublattice");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.4sl.toml";
    auto fl = FileToml(lfile);
    auto group = fl["Symmetries"].as<PermutationGroup>();
    auto action = GroupAction(group);
    auto action_sl = GroupActionSublattice<bit_t, 4>(group);
  }

  // Five sublattice
  {
    Log("GroupActionSublattice: 5 sublattice");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.10.heisenberg.5sl.toml";
    auto fl = FileToml(lfile);
    auto group = fl["Symmetries"].as<PermutationGroup>();
    auto action = GroupAction(group);
    auto action_sl = GroupActionSublattice<bit_t, 5>(group);
  }
}

TEST_CASE("GroupActionSublattice", "[symmetries]") {
  Log("Test GroupActionSublattice");
  Log("uint16_t");
  test_group_action_sublattice<uint16_t>();
  Log("uint32_t");
  test_group_action_sublattice<uint32_t>();
  Log("uint64_t");
  test_group_action_sublattice<uint64_t>();
  Log("Done");
}
