#include "../../catch.hpp"

#include <iostream>

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/symmetries/group_action/group_action.h>
#include <hydra/symmetries/group_action/group_action_lookup.h>

using namespace hydra;

template <class bit_t> void test_permutation_group_lookup(int n_sites) {
  using combinatorics::Subsets;

  // test cyclic group
  std::vector<Permutation> permutation_array;
  for (int sym = 0; sym < n_sites; ++sym) {

    std::vector<int> pv;
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      pv.push_back(newsite);
    }
    permutation_array.push_back(Permutation(pv));
  }
  auto perm_group = PermutationGroup(permutation_array);
  auto action = GroupAction(perm_group);
  auto lookup = GroupActionLookup<bit_t>(perm_group);

  // Check whether representative is smallest in orbit
  for (int sym = 0; sym < perm_group.size(); ++sym)
    for (bit_t state : Subsets<bit_t>(n_sites)) {

      auto s1 = action.apply(sym, state);
      auto s2 = lookup.apply(sym, state);
      CHECK(s1 == s2);
    }
}

TEST_CASE("GroupActionLookup", "[symmetries]") {
  hydra::Log.out("PermutationGroupLookup <-> GroupAction comparison");
  for (int n_sites = 1; n_sites < 6; ++n_sites) {
    test_permutation_group_lookup<uint16_t>(n_sites);
    test_permutation_group_lookup<uint32_t>(n_sites);
    test_permutation_group_lookup<uint64_t>(n_sites);
  }
  hydra::Log("done");
}
