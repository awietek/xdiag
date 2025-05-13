// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <iostream>

#include <xdiag/common.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/symmetries/operations/fermi_sign.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

template <class bit_t> void test_fermi_sign(int64_t nsites) {
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

  // auto sym_op = SpaceGroupOperator<bit_t>(nsites, permutation_array);

  auto fermi_work = symmetries::fermi_work(nsites);
  auto fermi_work_sort = symmetries::fermi_work_sort(nsites);
  for (auto state : Subsets<bit_t>(nsites)) {
    for (int64_t sym = 0; sym < nsites; ++sym) {

      // auto tstate = sym_op.apply(sym, state);
      // std::cout << "sym " << sym << "\n";
      // std::cout << bits_to_string(state, nsites) << "\n";
      // std::cout << bits_to_string(tstate, nsites) << "\n";

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
  xdiag::Log("Test fermi_sign");

  for (int64_t nsites = 1; nsites < 8; ++nsites) {
    test_fermi_sign<uint16_t>(nsites);
    test_fermi_sign<uint32_t>(nsites);
    test_fermi_sign<uint64_t>(nsites);
  }
  xdiag::Log("done");
}
