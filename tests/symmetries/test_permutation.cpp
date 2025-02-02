#include "../catch.hpp"

#include <iostream>
#include <random>

#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/utils/logger.hpp>

template <typename bit_t> void test_permutation_apply(int64_t nsites) {

  using namespace xdiag;
  std::random_device rd;
  std::mt19937 g(rd());
  for (int64_t i = 0; i < 20; ++i) {
    auto id = Permutation(nsites);
    auto a = id.array();
    std::shuffle(a.begin(), a.end(), g);
    auto p = Permutation(a);

    for (auto state : combinatorics::Subsets<bit_t>(nsites)) {
      bit_t tstate = 0;
      for (int64_t site = 0; site < nsites; ++site) {
        tstate |= ((state >> site) & 1) << p[site];
      }
      REQUIRE(tstate == p.apply(state));
    }
  }
}

TEST_CASE("permutation", "[symmetries]") {
  using namespace xdiag;
  Log("Test Permutation");

  std::random_device rd;
  std::mt19937 g(rd());
 
  // Test if identity is correct
  for (int64_t nsites = 1; nsites < 8; ++nsites) {
    std::vector<int64_t> pv(nsites);
    for (int64_t i = 0; i < nsites; ++i) {
      pv[i] = i;
    }
    auto p1 = Permutation(pv);
    auto p2 = Permutation(nsites);
    REQUIRE(p1 == p2);
  }

  for (int64_t nsites = 1; nsites < 8; ++nsites) {

    // Test identity multiplies
    for (int64_t i = 0; i < 5; ++i) {
      auto id = Permutation(nsites);
      auto a = id.array();
      std::shuffle(a.begin(), a.end(), g);
      auto p1 = Permutation(a);
      auto pi = Permutation(nsites);
      REQUIRE(p1 * pi == p1);
      REQUIRE(pi * p1 == p1);
    }

    // Test inverse
    for (int64_t i = 0; i < 20; ++i) {
      auto id = Permutation(nsites);
      auto a = id.array();
      std::shuffle(a.begin(), a.end(), g);
      auto p = Permutation(a);
      auto pinv = inverse(p);
      REQUIRE(p * pinv == id);
    }

    // Test applying permutations to state
    test_permutation_apply<uint16_t>(nsites);
    test_permutation_apply<uint32_t>(nsites);
    test_permutation_apply<uint64_t>(nsites);
  }

  Log("done");
}
