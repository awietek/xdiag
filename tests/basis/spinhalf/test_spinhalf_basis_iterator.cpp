#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"
#include <iostream>
#include <xdiag/blocks/spinhalf.hpp>

using namespace xdiag;

void test_spinhalf_basis_iterator(Spinhalf const &block) {

  int64_t n_sites = block.n_sites();
  auto pprev = ProductState(n_sites);
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
  for (int n_sites = 1; n_sites <= 6; ++n_sites) {
    auto block = Spinhalf(n_sites);
    test_spinhalf_basis_iterator(block);
  }

  Log("Spinhalf Basis Iterator Sz");
  for (int n_sites = 1; n_sites <= 6; ++n_sites) {
    for (int n_up = 0; n_up <= n_sites; ++n_up) {
      auto block = Spinhalf(n_sites, n_up);
      test_spinhalf_basis_iterator(block);
    }
  }

  Log("Spinhalf Basis Iterator No Sz Symmetric");
  for (int n_sites = 1; n_sites <= 6; ++n_sites) {
    auto irreps = testcases::electron::get_cyclic_group_irreps(n_sites);
    for (auto irrep : irreps) {
      auto block = Spinhalf(n_sites, irrep);
      test_spinhalf_basis_iterator(block);
    }
  }

  Log("Spinhalf Basis Iterator Sz Symmetric");
  for (int n_sites = 1; n_sites <= 6; ++n_sites) {
    auto irreps = testcases::electron::get_cyclic_group_irreps(n_sites);
    for (int n_up = 0; n_up <= n_sites; ++n_up) {
      for (auto irrep : irreps) {
        auto block = Spinhalf(n_sites, irrep);
        test_spinhalf_basis_iterator(block);
      }
    }
  }
}
