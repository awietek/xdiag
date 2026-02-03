// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/bitvector.hpp>
#include <xdiag/utils/logger.hpp>

#include <algorithm>
#include <limits>
#include <random>
#include <vector>

using namespace xdiag;
using namespace xdiag::bits;

template <typename bit_t> void test_bitvector_basic() {
  std::mt19937 rng(42);

  for (int nbits = 1; nbits <= std::numeric_limits<bit_t>::digits; ++nbits) {
    for (int size = 0; size < 100; ++size) {
      auto vec = BitVector<bit_t>(size, nbits);

      REQUIRE(vec.size() == size);
      REQUIRE(vec.nbits() == nbits);
      REQUIRE(vec.empty() == (size == 0));

      std::uniform_int_distribution<bit_t> dist(0, bitmask<bit_t>(nbits));
      for (int idx = 0; idx < size; ++idx) {
        for (int r = 0; r < 4; ++r) {
          bit_t bits = dist(rng);
          vec[idx] = bits;
          bit_t bits2 = vec[idx];
          REQUIRE(bits == bits2);
        }
      }
    }
  }
}

template <typename bit_t> void test_bitvector_element_access() {
  std::mt19937 rng(42);
  std::uniform_int_distribution<bit_t> dist(0, 255);

  // Test at() with bounds checking
  auto vec = BitVector<bit_t>(10, 8);

  for (int i = 0; i < 10; ++i) {
    bit_t val = dist(rng);
    vec.at(i) = val;
    REQUIRE(vec.at(i) == val);
  }

  // Test out of bounds throws
  REQUIRE_THROWS(vec.at(-1));
  REQUIRE_THROWS(vec.at(10));

  // Test front() and back()
  if (vec.size() > 0) {
    vec[0] = 42;
    vec[vec.size() - 1] = 99;
    REQUIRE(vec.front() == 42);
    REQUIRE(vec.back() == 99);
  }
}

template <typename bit_t> void test_bitvector_iterators() {
  std::mt19937 rng(42);

  auto vec = BitVector<bit_t>(20, 8);
  std::vector<bit_t> expected;

  std::uniform_int_distribution<bit_t> dist(0, 255);

  // Fill with random values
  for (int i = 0; i < 20; ++i) {
    bit_t val = dist(rng);
    vec[i] = val;
    expected.push_back(val);
  }

  // Test const iteration
  {
    std::vector<bit_t> result;
    for (auto it = vec.cbegin(); it != vec.cend(); ++it) {
      result.push_back(*it);
    }
    REQUIRE(result == expected);
  }

  // Test range-based for loop (const)
  {
    std::vector<bit_t> result;
    BitVector<bit_t> const &cvec = vec;
    for (auto val : cvec) {
      result.push_back(val);
    }
    REQUIRE(result == expected);
  }

  // Test non-const iteration
  {
    for (auto it = vec.begin(); it != vec.end(); ++it) {
      *it = (*it) ^ 1; // Flip lowest bit
    }
    for (size_t i = 0; i < expected.size(); ++i) {
      REQUIRE(vec[i] == (expected[i] ^ 1));
    }
  }

  // Test iterator arithmetic
  {
    auto it1 = vec.begin();
    auto it2 = it1 + 5;
    REQUIRE(it2 - it1 == 5);
    REQUIRE(*it2 == vec[5]);

    it2 += 3;
    REQUIRE(*it2 == vec[8]);

    it2 -= 2;
    REQUIRE(*it2 == vec[6]);
  }

  // Test STL algorithm (std::find)
  {
    bit_t search_val = vec[10];
    auto it = std::find(vec.begin(), vec.end(), search_val);
    REQUIRE(it != vec.end());
    REQUIRE(*it == search_val);
  }
}

template <typename bit_t> void test_bitvector_comparison() {
  auto vec1 = BitVector<bit_t>(10, 8);
  auto vec2 = BitVector<bit_t>(10, 8);
  auto vec3 = BitVector<bit_t>(10, 7); // different nbits

  for (int i = 0; i < 10; ++i) {
    vec1[i] = i;
    vec2[i] = i;
  }

  REQUIRE(vec1 == vec2);
  REQUIRE(!(vec1 != vec2));

  vec2[5] = 99;
  REQUIRE(vec1 != vec2);
  REQUIRE(!(vec1 == vec2));

  REQUIRE(vec1 != vec3); // different nbits
}

template <typename bit_t> void test_bitvector_max_nbits() {
  // Test with maximum nbits (equal to chunk size)
  constexpr int max_nbits = std::numeric_limits<bit_t>::digits;
  constexpr bit_t max_value = std::numeric_limits<bit_t>::max();

  auto vec = BitVector<bit_t>(100, max_nbits);

  REQUIRE(vec.size() == 100);
  REQUIRE(vec.nbits() == max_nbits);

  std::mt19937 rng(42);
  std::uniform_int_distribution<bit_t> dist(
      std::numeric_limits<bit_t>::min(),
      std::numeric_limits<bit_t>::max());

  // Test setting and getting full-width values
  for (int i = 0; i < 100; ++i) {
    bit_t val = dist(rng);
    vec[i] = val;
    REQUIRE(vec[i] == val);
  }

  // Test with max value
  vec[0] = max_value;
  REQUIRE(vec[0] == max_value);

  // Test with zero
  vec[1] = 0;
  REQUIRE(vec[1] == 0);

  // Test iteration with max nbits
  std::vector<bit_t> expected;
  for (int i = 0; i < 50; ++i) {
    bit_t val = dist(rng);
    vec[i] = val;
    expected.push_back(val);
  }

  int idx = 0;
  for (auto it = vec.begin(); it != vec.begin() + 50; ++it, ++idx) {
    REQUIRE(*it == expected[idx]);
  }

  // Test that nbits > max_nbits throws in constructor
  REQUIRE_THROWS(BitVector<bit_t>(10, max_nbits + 1));
}

TEST_CASE("bitvector", "[bits]") try {
  Log("Testing bitvector basic operations");

  SECTION("basic") {
    test_bitvector_basic<uint8_t>();
    test_bitvector_basic<uint16_t>();
    test_bitvector_basic<uint32_t>();
    test_bitvector_basic<uint64_t>();
  }

  SECTION("element_access") {
    test_bitvector_element_access<uint8_t>();
    test_bitvector_element_access<uint16_t>();
    test_bitvector_element_access<uint32_t>();
    test_bitvector_element_access<uint64_t>();
  }

  SECTION("iterators") {
    test_bitvector_iterators<uint8_t>();
    test_bitvector_iterators<uint16_t>();
    test_bitvector_iterators<uint32_t>();
    test_bitvector_iterators<uint64_t>();
  }

  SECTION("comparison") {
    test_bitvector_comparison<uint8_t>();
    test_bitvector_comparison<uint16_t>();
    test_bitvector_comparison<uint32_t>();
    test_bitvector_comparison<uint64_t>();
  }

  SECTION("max_nbits") {
    test_bitvector_max_nbits<uint8_t>();
    test_bitvector_max_nbits<uint16_t>();
    test_bitvector_max_nbits<uint32_t>();
    test_bitvector_max_nbits<uint64_t>();
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
