#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

template <typename bit_t> void test_permutation_apply(int n_sites) {

  using namespace hydra;

  for (int i = 0; i < 20; ++i) {
    auto id = identity_permutation(n_sites);
    auto p = shuffle(id);

    for (auto state : combinatorics::Subsets<bit_t>(n_sites)) {
      bit_t tstate = 0;
      for (int site = 0; site < n_sites; ++site) {
        tstate |= ((state >> site) & 1) << p[site];
      }
      REQUIRE(tstate == p.apply(state));
    }
  }
}

TEST_CASE("permutation", "[symmetries]") {
  using namespace hydra;
  Log("Test Permutation");

  // Test if identity is correct
  for (int n_sites = 1; n_sites < 8; ++n_sites) {
    std::vector<int> pv(n_sites);
    for (int i = 0; i < n_sites; ++i) {
      pv[i] = i;
    }
    auto p1 = Permutation(pv);
    auto p2 = identity_permutation(n_sites);
    REQUIRE(p1 == p2);
  }

  for (int n_sites = 1; n_sites < 8; ++n_sites) {

    // Test identity multiplies
    for (int i = 0; i < 5; ++i) {
      auto p1 = identity_permutation(n_sites).shuffle();
      auto pi = identity_permutation(n_sites);
      REQUIRE(p1 * pi == p1);
      REQUIRE(pi * p1 == p1);
    }

    // Test inverse
    for (int i = 0; i < 20; ++i) {
      auto id = identity_permutation(n_sites);
      auto p = shuffle(id);
      auto pinv = inverse(p);
      REQUIRE(p * pinv == id);
    }

    // Test applying permutations to state
    test_permutation_apply<uint16_t>(n_sites);
    test_permutation_apply<uint32_t>(n_sites);
    test_permutation_apply<uint64_t>(n_sites);
  }

  Log("done");
}
