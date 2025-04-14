// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"
#include <iostream>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

void test_spinhalf_basis_iterator(Spinhalf const &block) {

  int64_t nsites = block.nsites();
  auto pprev = ProductState(nsites);
  int64_t idx = 0;
  for (auto const &p : block) {
    REQUIRE(p != pprev);
    pprev = p;
    int64_t idx2 = block.index(p);
    // Log("{} {} {}", to_string(p), idx, idx2);
    REQUIRE(idx == idx2);
    ++idx;
  }
  REQUIRE(idx == block.dim());
}

TEST_CASE("spinhalf_basis_iterator", "[basis]") {

  Log("Spinhalf Basis Iterator No Sz");
  for (int nsites = 1; nsites <= 6; ++nsites) {
    auto block = Spinhalf(nsites);
    test_spinhalf_basis_iterator(block);
  }

  Log("Spinhalf Basis Iterator Sz");
  for (int nsites = 1; nsites <= 6; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      auto block = Spinhalf(nsites, nup);
      test_spinhalf_basis_iterator(block);
    }
  }

  Log("Spinhalf Basis Iterator No Sz Symmetric");
  for (int nsites = 1; nsites <= 6; ++nsites) {
    auto irreps = testcases::electron::get_cyclic_group_irreps(nsites);
    for (auto irrep : irreps) {
      auto block = Spinhalf(nsites, irrep);
      test_spinhalf_basis_iterator(block);
    }
  }

  Log("Spinhalf Basis Iterator Sz Symmetric");
  for (int nsites = 1; nsites <= 6; ++nsites) {
    auto irreps = testcases::electron::get_cyclic_group_irreps(nsites);
    for (int nup = 0; nup <= nsites; ++nup) {
      for (auto irrep : irreps) {
        auto block = Spinhalf(nsites, irrep);
        test_spinhalf_basis_iterator(block);
      }
    }
  }
}
