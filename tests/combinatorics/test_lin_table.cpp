// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <vector>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/utils/logger.hpp>

// Iterate LinTable and verify ordering, popcount, and size.
template <class bit_t> void test_lintable_iteration() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  for (int n = 0; n < 10; ++n) {
    for (int k = 0; k <= n; ++k) {
      LinTable<bit_t> lt(n, k);
      REQUIRE(lt.n() == n);
      REQUIRE(lt.k() == k);
      int64_t ctr = 0;
      bit_t current = 0;
      for (auto bits : lt) {
        if (ctr != 0)
          REQUIRE(bits > current);
        current = bits;
        ++ctr;
        REQUIRE(popcount(bits) == k);
      }
      REQUIRE(ctr == lt.size());
    }
  }
}

// index() must equal the sequential position in iteration.
template <class bit_t> void test_lintable_index() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 0; n < 10; ++n) {
    for (int k = 0; k <= n; ++k) {
      LinTable<bit_t> lt(n, k);
      int64_t idx = 0;
      for (auto bits : lt) {
        REQUIRE(lt.index(bits) == idx);
        ++idx;
      }
    }
  }
}

// operator[](i) matches sequential order; index() is its inverse.
template <class bit_t> void test_lintable_random_access() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 0; n < 10; ++n) {
    for (int k = 0; k <= n; ++k) {
      LinTable<bit_t> lt(n, k);

      std::vector<bit_t> elems;
      for (auto bits : lt)
        elems.push_back(bits);

      for (int64_t i = 0; i < lt.size(); ++i)
        REQUIRE(lt[i] == elems[i]);

      for (int64_t i = 0; i < lt.size(); ++i)
        REQUIRE(lt.index(elems[i]) == i);
    }
  }
}

// iterator operator+ and operator+=
template <class bit_t> void test_lintable_iterator_advance() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 1; n < 10; ++n) {
    for (int k = 0; k <= n; ++k) {
      LinTable<bit_t> lt(n, k);
      if (lt.size() == 0)
        continue;

      std::vector<bit_t> elems;
      for (auto bits : lt)
        elems.push_back(bits);

      for (int64_t i = 0; i < lt.size(); ++i)
        REQUIRE(*(lt.begin() + i) == elems[i]);

      auto it = lt.begin();
      for (int64_t i = 0; i < lt.size() - 1; ++i) {
        it += 1;
        REQUIRE(*it == elems[i + 1]);
      }

      if (lt.size() >= 3) {
        auto it2 = lt.begin();
        it2 += 2;
        REQUIRE(*it2 == elems[2]);
      }
    }
  }
}

// LinTable must produce the exact same sequence as Combinations.
template <class bit_t> void test_lintable_vs_combinations() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 0; n < 10; ++n) {
    for (int k = 0; k <= n; ++k) {
      LinTable<bit_t> lt(n, k);
      Combinations<bit_t> combs(n, k);

      REQUIRE(lt.n() == combs.n());
      REQUIRE(lt.k() == combs.k());
      REQUIRE(lt.size() == combs.size());

      auto it_lt = lt.begin();
      auto it_c = combs.begin();
      while (it_lt != lt.end()) {
        REQUIRE(*it_lt == *it_c);
        ++it_lt;
        ++it_c;
      }
      REQUIRE(it_c == combs.end());
    }
  }
}

TEST_CASE("LinTable", "[combinatorics]") {

  SECTION("iteration") {
    xdiag::Log("Testing LinTable - iteration");
    test_lintable_iteration<uint16_t>();
    test_lintable_iteration<uint32_t>();
    test_lintable_iteration<uint64_t>();
  }

  SECTION("index") {
    xdiag::Log("Testing LinTable - index");
    test_lintable_index<uint16_t>();
    test_lintable_index<uint32_t>();
    test_lintable_index<uint64_t>();
  }

  SECTION("random access and index round-trip") {
    xdiag::Log("Testing LinTable - random access and index round-trip");
    test_lintable_random_access<uint16_t>();
    test_lintable_random_access<uint32_t>();
    test_lintable_random_access<uint64_t>();
  }

  SECTION("iterator advance (+ and +=)") {
    xdiag::Log("Testing LinTable - iterator advance");
    test_lintable_iterator_advance<uint16_t>();
    test_lintable_iterator_advance<uint32_t>();
    test_lintable_iterator_advance<uint64_t>();
  }

  SECTION("vs Combinations") {
    xdiag::Log("Testing LinTable - vs Combinations");
    test_lintable_vs_combinations<uint16_t>();
    test_lintable_vs_combinations<uint32_t>();
    test_lintable_vs_combinations<uint64_t>();
  }
}
