#include "../catch.hpp"

#include <iostream>

#include <hydra/common.h>
#include <hydra/symmetries/permutation.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/symmetries/operations/fermi_sign.h>

using namespace hydra;

template <class bit_t> void test_fermi_sign(int64_t n_sites) {
  using combinatorics::Subsets;

  // test cyclic group
  std::vector<Permutation> permutation_array;
  for (int64_t sym = 0; sym < n_sites; ++sym) {

    std::vector<int64_t> pv;
    for (int64_t site = 0; site < n_sites; ++site) {
      int64_t newsite = (site + sym) % n_sites;
      pv.push_back(newsite);
    }
    permutation_array.push_back(Permutation(pv));
  }

  // auto sym_op = SpaceGroupOperator<bit_t>(n_sites, permutation_array);

  auto fermi_work = symmetries::fermi_work(n_sites);
  auto fermi_work_sort = symmetries::fermi_work_sort(n_sites);
  for (auto state : Subsets<bit_t>(n_sites)) {
    for (int64_t sym = 0; sym < n_sites; ++sym) {

      // auto tstate = sym_op.apply(sym, state);
      // std::cout << "sym " << sym << "\n";
      // std::cout << bits_to_string(state, n_sites) << "\n";
      // std::cout << bits_to_string(tstate, n_sites) << "\n";

      auto fermi1 = symmetries::fermi_sign_of_permutation(
          state, permutation_array[sym], fermi_work);
      auto fermi2 = symmetries::fermi_sign_of_permutation_sort(
          state, permutation_array[sym], fermi_work_sort);

      // std::cout << fermi1 << " " << fermi2 << "\n";
      // std::cout << "\n";
      REQUIRE(fermi1 == fermi2);
    }
    // std::cout << "-----\n";
  }
}

TEST_CASE("fermi_sign", "[symmetries]") {
  hydra::Log("Test fermi_sign");

  for (int64_t n_sites = 1; n_sites < 8; ++n_sites) {
    test_fermi_sign<uint16_t>(n_sites);
    test_fermi_sign<uint32_t>(n_sites);
    test_fermi_sign<uint64_t>(n_sites);
  }
  hydra::Log("done");
}
