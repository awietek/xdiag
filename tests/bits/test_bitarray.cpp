// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/utils/logger.hpp>

#include <random>

using namespace xdiag;
using namespace xdiag::bits;

template <typename bit_t, int nbits>
void test_bitarray_basic() {
  BitArray<bit_t, nbits> arr;

  // Test default initialization (all zeros)
  for (int64_t i = 0; i < BitArray<bit_t, nbits>::maximum_size; ++i) {
    REQUIRE(arr.get(i) == 0);
  }

  // Test setting and getting values
  int64_t max_value = (1LL << nbits) - 1;
  for (int64_t i = 0; i < std::min(int64_t(10), BitArray<bit_t, nbits>::maximum_size); ++i) {
    int64_t value = (i * 3) & max_value; // Some pattern within range
    arr.set(i, value);
    REQUIRE(arr.get(i) == value);
  }
}

template <typename bit_t, int nbits>
void test_bitarray_boundaries() {
  BitArray<bit_t, nbits> arr;
  int64_t max_value = (1LL << nbits) - 1;

  // Test minimum and maximum values
  arr.set(0, 0);
  REQUIRE(arr.get(0) == 0);

  arr.set(0, max_value);
  REQUIRE(arr.get(0) == max_value);

  // Test that values wrap/truncate correctly
  arr.set(0, max_value + 1);
  REQUIRE(arr.get(0) == 0); // Should wrap to 0

  arr.set(0, -1);
  REQUIRE(arr.get(0) == max_value); // -1 has all bits set
}

template <typename bit_t, int nbits>
void test_bitarray_independence() {
  BitArray<bit_t, nbits> arr;
  int64_t max_value = (1LL << nbits) - 1;
  int64_t size = std::min(int64_t(8), BitArray<bit_t, nbits>::maximum_size);

  // Set different values at different positions
  for (int64_t i = 0; i < size; ++i) {
    arr.set(i, i & max_value);
  }

  // Verify each position retained its value
  for (int64_t i = 0; i < size; ++i) {
    REQUIRE(arr.get(i) == (i & max_value));
  }

  // Modify one position and verify others unchanged
  if (size > 2) {
    arr.set(1, max_value);
    REQUIRE(arr.get(0) == 0);
    REQUIRE(arr.get(1) == max_value);
    REQUIRE(arr.get(2) == 2);
  }
}

template <typename bit_t, int nbits>
void test_bitarray_random() {
  std::mt19937 rng(42 + nbits);
  std::uniform_int_distribution<int64_t> dist(0, (1LL << nbits) - 1);

  BitArray<bit_t, nbits> arr;
  int64_t size = std::min(int64_t(20), BitArray<bit_t, nbits>::maximum_size);

  // Set random values
  std::vector<int64_t> expected(size);
  for (int64_t i = 0; i < size; ++i) {
    expected[i] = dist(rng);
    arr.set(i, expected[i]);
  }

  // Verify all values
  for (int64_t i = 0; i < size; ++i) {
    REQUIRE(arr.get(i) == expected[i]);
  }
}

template <typename bit_t, int nbits>
void test_bitarray_comparison() {
  BitArray<bit_t, nbits> arr1, arr2;

  // Test equality of default-constructed arrays
  REQUIRE(arr1 == arr2);
  REQUIRE(!(arr1 != arr2));

  // Modify one array
  arr1.set(0, 1);
  REQUIRE(arr1 != arr2);
  REQUIRE(!(arr1 == arr2));

  // Make them equal again
  arr2.set(0, 1);
  REQUIRE(arr1 == arr2);
  REQUIRE(!(arr1 != arr2));
}

TEST_CASE("bitarray", "[bits]") {

  SECTION("basic operations") {
    Log("Testing BitArray - basic operations");
    // Native integers
    test_bitarray_basic<uint16_t, 1>();
    test_bitarray_basic<uint16_t, 2>();
    test_bitarray_basic<uint16_t, 4>();
    test_bitarray_basic<uint32_t, 3>();
    test_bitarray_basic<uint32_t, 5>();
    test_bitarray_basic<uint64_t, 6>();
    test_bitarray_basic<uint64_t, 8>();

    // Bitsets
    test_bitarray_basic<Bitset<uint64_t, 1>, 2>();
    test_bitarray_basic<Bitset<uint64_t, 2>, 4>();
    test_bitarray_basic<Bitset<uint64_t, 4>, 6>();
  }

  SECTION("boundary values") {
    Log("Testing BitArray - boundary values");
    test_bitarray_boundaries<uint16_t, 1>();
    test_bitarray_boundaries<uint16_t, 4>();
    test_bitarray_boundaries<uint32_t, 3>();
    test_bitarray_boundaries<uint64_t, 7>();

    test_bitarray_boundaries<Bitset<uint64_t, 1>, 2>();
    test_bitarray_boundaries<Bitset<uint64_t, 2>, 5>();
  }

  SECTION("element independence") {
    Log("Testing BitArray - element independence");
    test_bitarray_independence<uint16_t, 2>();
    test_bitarray_independence<uint32_t, 3>();
    test_bitarray_independence<uint64_t, 4>();

    test_bitarray_independence<Bitset<uint64_t, 1>, 3>();
    test_bitarray_independence<Bitset<uint64_t, 2>, 5>();
  }

  SECTION("random values") {
    Log("Testing BitArray - random values");
    test_bitarray_random<uint16_t, 1>();
    test_bitarray_random<uint16_t, 3>();
    test_bitarray_random<uint32_t, 5>();
    test_bitarray_random<uint64_t, 7>();

    test_bitarray_random<Bitset<uint64_t, 1>, 2>();
    test_bitarray_random<Bitset<uint64_t, 2>, 4>();
    test_bitarray_random<Bitset<uint64_t, 4>, 6>();
  }

  SECTION("comparison operators") {
    Log("Testing BitArray - comparison operators");
    test_bitarray_comparison<uint16_t, 2>();
    test_bitarray_comparison<uint32_t, 4>();
    test_bitarray_comparison<uint64_t, 6>();

    test_bitarray_comparison<Bitset<uint64_t, 1>, 3>();
    test_bitarray_comparison<Bitset<uint64_t, 2>, 5>();
  }
}
