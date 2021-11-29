#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

using namespace hydra;

template <class bit_t> void test_permutation_group_lookup(int n_sites) {
  using combinatorics::Subsets;

  // test cyclic group
  std::vector<int> permutation_array;
  for (int sym = 0; sym < n_sites; ++sym) {
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation_array.push_back(newsite);
    }
  }

  auto perm_group = PermutationGroup(n_sites, n_sites, permutation_array);
  auto action = PermutationGroupAction(perm_group);
  auto lookup = PermutationGroupLookup<bit_t>(perm_group);

  // Check whether representative is smallest in orbit
  for (int sym =0; sym< perm_group.size(); ++sym)
  for (bit_t state : Subsets<bit_t>(n_sites)) {

    auto s1 = action.apply(sym, state);
    auto s2 = lookup.apply(sym, state);
    CHECK(s1 == s2);
  }
}

TEST_CASE("PermutationGroupLookup", "[symmetries]") {
  lila::Log.out("PermutationGroupLookup <-> PermutationGroupAction comparison");
  for (int n_sites = 1; n_sites < 6; ++n_sites) {
    test_permutation_group_lookup<uint16_t>(n_sites);
    test_permutation_group_lookup<uint32_t>(n_sites);
    test_permutation_group_lookup<uint64_t>(n_sites);
  }
  lila::Log("done");
}
