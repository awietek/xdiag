#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

using namespace hydra;

template <class bit_t> void test_spacegroup_operator(int n_sites) {

  // test cyclic group
  std::vector<int> permutation_array;
  for (int sym = 0; sym < n_sites; ++sym) {
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation_array.push_back(newsite);
    }
  }

  auto sym_op = SpaceGroupOperator<bit_t>(n_sites, permutation_array);

  // // Print the action on the states
  // for (auto state : Subsets<bit_t>(n_sites)) {
  //   std::cout << bits_to_string(state, n_sites) << "\n";

  //   for (int sym = 0; sym < n_sites; ++sym) {
  //     auto tstate = sym_op.apply(sym, state);
  //     std::cout << sym << " " << bits_to_string(tstate, n_sites) << "\n";
  //   }
  //   std::cout << "\n";
  // }

  // Check whether cyclic group property holds for every state
  for (auto state : Subsets<bit_t>(n_sites)) {
    for (int sym1 = 0; sym1 < n_sites; ++sym1)
      for (int sym2 = 0; sym2 < n_sites; ++sym2) {
        auto tstate1 = sym_op.apply(sym1, state);
        auto tstate2 = sym_op.apply(sym2, tstate1);
        auto tstate12 = sym_op.apply((sym1 + sym2) % n_sites, state);
        REQUIRE(tstate12 == tstate2);
      }
  }

  // Check whether representative is smallest in orbit
  for (auto state : Subsets<bit_t>(n_sites)) {

    auto rep = sym_op.representative(state);

    // std::cout << "s " << bits_to_string(state, n_sites) << "\n";
    // std::cout << "r " << bits_to_string(rep, n_sites) << "\n\n";

    std::vector<bit_t> orbit;
    for (int sym = 0; sym < n_sites; ++sym) {
      orbit.push_back(sym_op.apply(sym, state));
    }
    auto min = std::min_element(orbit.begin(), orbit.end());
    REQUIRE(*min == rep);

    auto [rep2, sym] = sym_op.representative_index(state);
    REQUIRE(rep2 == rep);
    REQUIRE(rep == sym_op.apply(sym, state));

    auto [rep3, nsym, sym_ptr] = sym_op.representative_indices(state);
    REQUIRE(rep3 == rep);

    // std::cout << "s " << bits_to_string(state, n_sites) << "\n";
    // std::cout << "r " << bits_to_string(rep, n_sites) << "\n";
    for (int i=0; i<nsym; ++i){
      auto tstate = sym_op.apply(sym_ptr[i], state);
      REQUIRE(rep == tstate);
      // std::cout << "t " << bits_to_string(tstate, n_sites) << " " << sym_ptr[i] << "\n";
    }
    // std::cout << "\n";
  }
}

TEST_CASE("spacegroup_operator", "[symmetries]") {
  
  for (int n_sites = 1; n_sites < 6; ++n_sites) {
    test_spacegroup_operator<hydra::uint16>(n_sites);
    test_spacegroup_operator<hydra::uint32>(n_sites);
    test_spacegroup_operator<hydra::uint64>(n_sites);
  }
}
