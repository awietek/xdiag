// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>
#include <xdiag/utils/logger.hpp>
using namespace xdiag;

template <class bit_t> void test_permutation_group_lookup(int64_t nsites) {
  using combinatorics::Subsets;

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
  auto perm_group = PermutationGroup(permutation_array);
  auto action = GroupAction(perm_group);
  auto lookup = GroupActionLookup<bit_t>(perm_group);

  // Check whether representative is smallest in orbit
  for (int64_t sym = 0; sym < perm_group.size(); ++sym)
    for (bit_t state : Subsets<bit_t>(nsites)) {

      auto s1 = action.apply(sym, state);
      auto s2 = lookup.apply(sym, state);
      CHECK(s1 == s2);
    }
}

TEST_CASE("GroupActionLookup", "[symmetries]") {
  xdiag::Log("PermutationGroupLookup <-> GroupAction comparison");
  for (int64_t nsites = 1; nsites < 6; ++nsites) {
    test_permutation_group_lookup<uint16_t>(nsites);
    test_permutation_group_lookup<uint32_t>(nsites);
    test_permutation_group_lookup<uint64_t>(nsites);
  }
  xdiag::Log("done");
}
