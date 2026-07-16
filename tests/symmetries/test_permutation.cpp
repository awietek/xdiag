// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <iostream>
#include <random>

#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/utils/logger.hpp>

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

    // Test inverse (both sides)
    for (int64_t i = 0; i < 20; ++i) {
      auto id = Permutation(nsites);
      auto a = id.array();
      std::shuffle(a.begin(), a.end(), g);
      auto p = Permutation(a);
      auto pinv = inv(p);
      REQUIRE(p * pinv == id);
      REQUIRE(pinv * p == id);
    }

    // Test associativity
    for (int64_t i = 0; i < 5; ++i) {
      auto id = Permutation(nsites);
      auto a1 = id.array(); std::shuffle(a1.begin(), a1.end(), g);
      auto a2 = id.array(); std::shuffle(a2.begin(), a2.end(), g);
      auto a3 = id.array(); std::shuffle(a3.begin(), a3.end(), g);
      auto p1 = Permutation(a1), p2 = Permutation(a2), p3 = Permutation(a3);
      REQUIRE((p1 * p2) * p3 == p1 * (p2 * p3));
    }
  }

  // Test pow: cyclic shift of order n satisfies p^n == id
  for (int64_t n = 1; n < 8; ++n) {
    std::vector<int64_t> pv(n);
    for (int64_t i = 0; i < n; ++i)
      pv[i] = (i + 1) % n;
    auto p = Permutation(pv);
    auto id = Permutation(n);
    REQUIRE(pow(p, 0) == id);
    REQUIRE(pow(p, n) == id);
    REQUIRE(pow(p, -1) == inv(p));
    REQUIRE(pow(p, 1) == p);
  }

  // Test construction from pointer
  for (int64_t nsites = 1; nsites < 8; ++nsites) {
    auto id = Permutation(nsites);
    auto a = id.array();
    std::shuffle(a.begin(), a.end(), g);
    auto p1 = Permutation(a);
    auto p2 = Permutation(a.data(), nsites);
    REQUIRE(p1 == p2);
  }

  Log("done");
}
