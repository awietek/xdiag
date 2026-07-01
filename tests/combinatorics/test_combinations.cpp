// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <vector>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/enumerate_combinations.hpp>
#include <xdiag/utils/logger.hpp>

template <typename bit_t> void test_combinations() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  for (int n = 0; n < 7; ++n) {
    for (int k = 0; k <= n; ++k) {
      Combinations<bit_t> combs(n, k);
      REQUIRE(n == combs.n());
      REQUIRE(k == combs.k());
      int64_t ctr = 0;
      bit_t current = 0;
      for (auto comb : combs) {

        if (ctr != 0)
          REQUIRE(comb > current);
        current = comb;
        ++ctr;
        REQUIRE(popcount(comb) == k);
        REQUIRE(comb < ((bit_t)1 << n));
      }
      REQUIRE(ctr == combs.size());
    }
  }
}

// Test Bitset combinations against native uint64_t reference
template <typename chunk_t, int64_t nchunks> void test_combinations_bitset() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  // Calculate maximum n based on Bitset capacity
  constexpr int64_t chunk_bits = std::numeric_limits<chunk_t>::digits;
  constexpr int64_t max_bits = (nchunks == 0) ? 64 : (nchunks * chunk_bits);
  constexpr int max_n = std::min(int64_t(16), max_bits - 1);

  // Test various (n, k) pairs
  for (int n = 0; n <= max_n; ++n) {
    for (int k = 0; k <= n; ++k) {
      Combinations<Bitset<chunk_t, nchunks>> combs_bitset(n, k);
      Combinations<uint64_t> combs_uint64(n, k);

      REQUIRE(combs_bitset.n() == combs_uint64.n());
      REQUIRE(combs_bitset.k() == combs_uint64.k());
      REQUIRE(combs_bitset.size() == combs_uint64.size());

      // Compare each generated pattern
      auto it_bitset = combs_bitset.begin();
      auto it_uint64 = combs_uint64.begin();

      int64_t count = 0;
      while (it_bitset != combs_bitset.end() &&
             it_uint64 != combs_uint64.end()) {
        auto pattern_bitset = *it_bitset;
        auto pattern_uint64 = *it_uint64;

        // Convert Bitset to uint64_t for comparison
        uint64_t bitset_as_uint64 = to_uint64(pattern_bitset);
        REQUIRE(bitset_as_uint64 == pattern_uint64);

        // Verify bit count
        REQUIRE(pattern_bitset.count() == k);
        REQUIRE(popcount(pattern_uint64) == k);

        ++it_bitset;
        ++it_uint64;
        ++count;
      }

      REQUIRE(count == combs_bitset.size());
      REQUIRE(it_bitset == combs_bitset.end());
      REQUIRE(it_uint64 == combs_uint64.end());
    }
  }
}

// Test Bitset with larger n values
template <typename chunk_t, int64_t nchunks>
void test_combinations_bitset_large() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  // Test cases that exceed single uint64_t capacity (n > 64)
  if constexpr (nchunks >= 2) {
    for (int n : {65, 70, 80, 100}) {
      for (int k : {0, 1, 2, n - 1, n}) {
        if (k > n)
          continue;

        Combinations<Bitset<chunk_t, nchunks>> combs(n, k);

        REQUIRE(combs.n() == n);
        REQUIRE(combs.k() == k);

        int64_t count = 0;
        Bitset<chunk_t, nchunks> prev(n);

        for (auto pattern : combs) {
          // Verify bit count
          REQUIRE(pattern.count() == k);

          // Verify patterns are increasing
          if (count > 0) {
            REQUIRE(pattern > prev);
          }

          prev = pattern;
          ++count;
        }

        REQUIRE(count == combs.size());
      }
    }
  }
}

template <typename chunk_t, int64_t nchunks>
void test_combinations_random_access_bitset() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  for (int n = 0; n < 7; ++n) {
    for (int k = 0; k <= n; ++k) {
      Combinations<Bitset<chunk_t, nchunks>> combs(n, k);

      std::vector<Bitset<chunk_t, nchunks>> elems;
      for (auto c : combs)
        elems.push_back(c);

      for (int64_t i = 0; i < combs.size(); ++i)
        REQUIRE(combs[i] == elems[i]);

      for (int64_t i = 0; i < combs.size(); ++i)
        REQUIRE(combs.index(elems[i]) == i);
    }
  }
}

template <typename chunk_t, int64_t nchunks>
void test_combinations_iterator_advance_bitset() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  for (int n = 1; n < 7; ++n) {
    for (int k = 0; k <= n; ++k) {
      Combinations<Bitset<chunk_t, nchunks>> combs(n, k);
      if (combs.size() == 0)
        continue;

      std::vector<Bitset<chunk_t, nchunks>> elems;
      for (auto c : combs)
        elems.push_back(c);

      for (int64_t i = 0; i < combs.size(); ++i)
        REQUIRE(*(combs.begin() + i) == elems[i]);

      auto it = combs.begin();
      for (int64_t i = 0; i < combs.size() - 1; ++i) {
        it += 1;
        REQUIRE(*it == elems[i + 1]);
      }

      if (combs.size() >= 3) {
        auto it2 = combs.begin();
        it2 += 2;
        REQUIRE(*it2 == elems[2]);
      }
    }
  }
}

template <typename bit_t> void test_combinations_random_access() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 0; n < 7; ++n) {
    for (int k = 0; k <= n; ++k) {
      Combinations<bit_t> combs(n, k);

      // Collect elements sequentially
      std::vector<bit_t> elems;
      for (auto c : combs)
        elems.push_back(c);

      // operator[]: combs[i] must match sequential order
      for (int64_t i = 0; i < combs.size(); ++i)
        REQUIRE(combs[i] == elems[i]);

      // index: round-trip index(combs[i]) == i
      for (int64_t i = 0; i < combs.size(); ++i)
        REQUIRE(combs.index(elems[i]) == i);
    }
  }
}

template <typename bit_t> void test_combinations_iterator_advance() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 1; n < 7; ++n) {
    for (int k = 0; k <= n; ++k) {
      Combinations<bit_t> combs(n, k);
      if (combs.size() == 0)
        continue;

      std::vector<bit_t> elems;
      for (auto c : combs)
        elems.push_back(c);

      // operator+: begin() + i matches elems[i]
      for (int64_t i = 0; i < combs.size(); ++i)
        REQUIRE(*(combs.begin() + i) == elems[i]);

      // operator+=: step-by-step advance
      auto it = combs.begin();
      for (int64_t i = 0; i < combs.size() - 1; ++i) {
        it += 1;
        REQUIRE(*it == elems[i + 1]);
      }

      // operator+=: larger jump
      if (combs.size() >= 3) {
        auto it2 = combs.begin();
        it2 += 2;
        REQUIRE(*it2 == elems[2]);
      }
    }
  }
}

// Verify next_combination(v, n) produces the same sequence as Combinations,
// and for integral bit_t also agrees with the one-arg next_combination(v).
template <typename bit_t> void test_next_combination_overload() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  for (int n = 1; n < 8; ++n) {
    for (int k = 1; k <= n; ++k) {
      // Collect reference sequence from Combinations
      std::vector<bit_t> elems;
      for (auto c : Combinations<bit_t>(n, k))
        elems.push_back(c);

      // Walk using next_combination(v, n) and compare
      bit_t v = elems[0];
      for (int64_t i = 0; i < (int64_t)elems.size() - 1; ++i) {
        bit_t got = next_combination(v, n);
        REQUIRE(got == elems[i + 1]);
        v = got;
      }

      // For integral types: two-arg must agree with one-arg on all elements
      // except the last (calling next_combination on the last pattern is UB
      // for the algorithm so we skip it)
      if constexpr (std::is_integral_v<bit_t>) {
        for (int64_t i = 0; i < (int64_t)elems.size() - 1; ++i)
          REQUIRE(next_combination(elems[i], n) == next_combination(elems[i]));
      }
    }
  }
}

template <typename chunk_t, int64_t nchunks>
void test_next_combination_overload_bitset() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  for (int n = 1; n < 8; ++n) {
    for (int k = 1; k <= n; ++k) {
      std::vector<Bitset<chunk_t, nchunks>> elems;
      for (auto c : Combinations<Bitset<chunk_t, nchunks>>(n, k))
        elems.push_back(c);

      Bitset<chunk_t, nchunks> v = elems[0];
      for (int64_t i = 0; i < (int64_t)elems.size() - 1; ++i) {
        Bitset<chunk_t, nchunks> got = next_combination(v, (int64_t)n);
        REQUIRE(got == elems[i + 1]);
        v = got;
      }
    }
  }
}

TEST_CASE("Combinations", "[combinatorics]") {
  using namespace xdiag::bits;

  SECTION("native integers") {
    xdiag::Log("Testing Combinations - native integers");
    test_combinations<uint16_t>();
    test_combinations<uint32_t>();
    test_combinations<uint64_t>();
  }

  SECTION("random access and index") {
    xdiag::Log("Testing Combinations - random access and index");
    test_combinations_random_access<uint16_t>();
    test_combinations_random_access<uint32_t>();
    test_combinations_random_access<uint64_t>();
  }

  SECTION("iterator advance (+ and +=)") {
    xdiag::Log("Testing Combinations - iterator advance");
    test_combinations_iterator_advance<uint16_t>();
    test_combinations_iterator_advance<uint32_t>();
    test_combinations_iterator_advance<uint64_t>();
  }

  SECTION("next_combination(v, n) - integers") {
    xdiag::Log("Testing Combinations - next_combination (integers)");
    test_next_combination_overload<uint16_t>();
    test_next_combination_overload<uint32_t>();
    test_next_combination_overload<uint64_t>();
  }

  SECTION("next_combination(v, n) - bitset") {
    xdiag::Log("Testing Combinations - next_combination (bitset)");
    test_next_combination_overload_bitset<uint8_t, 0>();
    test_next_combination_overload_bitset<uint8_t, 1>();
    test_next_combination_overload_bitset<uint64_t, 1>();
    test_next_combination_overload_bitset<uint64_t, 2>();
  }

  SECTION("random access and index (bitset)") {
    xdiag::Log("Testing Combinations - random access and index (bitset)");
    test_combinations_random_access_bitset<uint8_t, 0>();
    test_combinations_random_access_bitset<uint8_t, 1>();
    test_combinations_random_access_bitset<uint64_t, 1>();
    test_combinations_random_access_bitset<uint64_t, 2>();
  }

  SECTION("iterator advance (+ and +=) (bitset)") {
    xdiag::Log("Testing Combinations - iterator advance (bitset)");
    test_combinations_iterator_advance_bitset<uint8_t, 0>();
    test_combinations_iterator_advance_bitset<uint8_t, 1>();
    test_combinations_iterator_advance_bitset<uint64_t, 1>();
    test_combinations_iterator_advance_bitset<uint64_t, 2>();
  }

  SECTION("bitset vs uint64_t") {
    xdiag::Log("Testing Bitset<uint8_t, 0> combinations");
    test_combinations_bitset<uint8_t, 0>();

    xdiag::Log("Testing Bitset<uint8_t, 1> combinations");
    test_combinations_bitset<uint8_t, 1>();

    xdiag::Log("Testing Bitset<uint16_t, 1> combinations");
    test_combinations_bitset<uint16_t, 1>();

    xdiag::Log("Testing Bitset<uint32_t, 1> combinations");
    test_combinations_bitset<uint32_t, 1>();

    xdiag::Log("Testing Bitset<uint64_t, 1> combinations");
    test_combinations_bitset<uint64_t, 1>();

    xdiag::Log("Testing Bitset<uint64_t, 2> combinations");
    test_combinations_bitset<uint64_t, 2>();
  }

  SECTION("bitset large n") {
    xdiag::Log("Testing Bitset<uint64_t, 2> large n");
    test_combinations_bitset_large<uint64_t, 2>();

    xdiag::Log("Testing Bitset<uint64_t, 4> large n");
    test_combinations_bitset_large<uint64_t, 4>();

    xdiag::Log("Testing Bitset<uint64_t, 8> large n");
    test_combinations_bitset_large<uint64_t, 8>();
  }
}
