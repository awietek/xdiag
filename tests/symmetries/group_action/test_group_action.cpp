#include "../../catch.hpp"

#include <iostream>

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>

using namespace xdiag;

template <class bit_t> void test_permutation_group_action(int64_t nsites) {
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
  auto sym_op = GroupAction(perm_group);

  // // Print the action on the states
  // for (auto state : Subsets<bit_t>(nsites)) {
  //   std::cout << bits_to_string(state, nsites) << "\n";

  //   for (int64_t sym = 0; sym < nsites; ++sym) {
  //     auto tstate = sym_op.apply(sym, state);
  //     std::cout << sym << " " << bits_to_string(tstate, nsites) << "\n";
  //   }
  //   std::cout << "\n";
  // }

  // Check whether cyclic group property holds for every state
  for (auto state : Subsets<bit_t>(nsites)) {
    for (int64_t sym1 = 0; sym1 < nsites; ++sym1)
      for (int64_t sym2 = 0; sym2 < nsites; ++sym2) {
        auto tstate1 = sym_op.apply(sym1, state);
        auto tstate2 = sym_op.apply(sym2, tstate1);
        auto tstate12 = sym_op.apply((sym1 + sym2) % nsites, state);
        REQUIRE(tstate12 == tstate2);
      }
  }

  // Check whether representative is smallest in orbit
  for (auto state : Subsets<bit_t>(nsites)) {

    auto rep = sym_op.representative(state);

    // std::cout << "s " << bits_to_string(state, nsites) << "\n";
    // std::cout << "r " << bits_to_string(rep, nsites) << "\n\n";

    std::vector<bit_t> orbit;
    for (int64_t sym = 0; sym < nsites; ++sym) {
      orbit.push_back(sym_op.apply(sym, state));
    }
    auto min = std::min_element(orbit.begin(), orbit.end());
    REQUIRE(*min == rep);

    auto [rep2, sym] = sym_op.representative_sym(state);
    REQUIRE(rep2 == rep);
    REQUIRE(rep == sym_op.apply(sym, state));

    auto [rep3, syms] = sym_op.representative_syms(state);
    REQUIRE(rep3 == rep);

    for (int64_t sym : syms) {
      auto tstate = sym_op.apply(sym, state);
      REQUIRE(rep == tstate);
    }
  }
}

TEST_CASE("GroupAction", "[symmetries]") {
  xdiag::Log("Test GroupAction");
  for (int64_t nsites = 1; nsites < 6; ++nsites) {
    test_permutation_group_action<uint16_t>(nsites);
    test_permutation_group_action<uint32_t>(nsites);
    test_permutation_group_action<uint64_t>(nsites);
  }
  xdiag::Log("done");
}
