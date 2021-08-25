#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

using namespace hydra;

template <class bit_t> void test_fermi_sign(int n_sites) {

  // test cyclic group
  std::vector<int> permutation_array;
  for (int sym = 0; sym < n_sites; ++sym) {
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation_array.push_back(newsite);
    }
  }

  // auto sym_op = SpaceGroupOperator<bit_t>(n_sites, permutation_array);

  std::vector<int> work(2 * n_sites);
  for (auto state : Subsets<bit_t>(n_sites)) {
    for (int sym = 0; sym < n_sites; ++sym) {

      // auto tstate = sym_op.apply(sym, state);
      // std::cout << "sym " << sym << "\n";
      // std::cout << bits_to_string(state, n_sites) << "\n";
      // std::cout << bits_to_string(tstate, n_sites) << "\n";

      auto fermi1 = hydra::utils::fermi_sign_of_permutation(
          state, permutation_array.data() + sym * n_sites, work.data());
      auto fermi2 = hydra::utils::fermi_sign_of_permutation_sort(
          state, permutation_array.data() + sym * n_sites, work.data());

      // std::cout << fermi1 << " " << fermi2 << "\n";
      // std::cout << "\n";
      REQUIRE(fermi1 == fermi2);
    }
    // std::cout << "-----\n";
  }
}

TEST_CASE("fermi_sign", "[symmetries]") {
  for (int n_sites = 1; n_sites < 8; ++n_sites) {
    test_fermi_sign<hydra::uint16>(n_sites);
    test_fermi_sign<hydra::uint32>(n_sites);
    test_fermi_sign<hydra::uint64>(n_sites);
  }
}
