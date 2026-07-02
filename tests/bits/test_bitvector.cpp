// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/bitvector.hpp>
#include <xdiag/utils/logger.hpp>

#include <algorithm>
#include <limits>
#include <random>
#include <vector>

using namespace xdiag;
using namespace xdiag::bits;

// --- Integral value_t (uint16_t, uint32_t, uint64_t) ---

template <typename value_t> void test_bitvector_integral() {
  std::mt19937 rng(42);

  for (int64_t nbits = 1; nbits <= std::numeric_limits<value_t>::digits;
       ++nbits) {
    value_t mask = (nbits == std::numeric_limits<value_t>::digits)
                       ? std::numeric_limits<value_t>::max()
                       : bitmask<value_t>(nbits);
    auto vec = BitVector<value_t>(20, nbits);
    REQUIRE(vec.size() == 20);
    REQUIRE(vec.nbits() == nbits);
    REQUIRE(!vec.empty());

    std::vector<value_t> expected;
    for (int i = 0; i < 20; ++i) {
      value_t val =
          std::uniform_int_distribution<value_t>(0, mask)(rng);
      vec[i] = val;
      expected.push_back(val);
    }
    for (int i = 0; i < 20; ++i)
      REQUIRE(vec[i] == expected[i]);

    // const iteration
    int idx = 0;
    for (auto v : (BitVector<value_t> const &)vec)
      REQUIRE(v == expected[idx++]);
  }

  // empty vector
  REQUIRE(BitVector<value_t>(0, 4).empty());

  // at() bounds checking
  {
    auto vec = BitVector<value_t>(5, 8);
    // REQUIRE_THROWS(vec.at(-1));
    // REQUIRE_THROWS(vec.at(5));

    // front/back
    vec[0] = 1;
    vec[4] = 7;
    REQUIRE(vec.front() == 1);
    REQUIRE(vec.back() == 7);
  }

  // operator== / !=
  {
    auto v1 = BitVector<value_t>(5, 8);
    auto v2 = BitVector<value_t>(5, 8);
    for (int i = 0; i < 5; ++i) {
      v1[i] = i;
      v2[i] = i;
    }
    REQUIRE(v1 == v2);
    v2[2] = 99;
    REQUIRE(v1 != v2);
    REQUIRE(BitVector<value_t>(5, 4) != BitVector<value_t>(5, 5));
  }
}

// --- BitsetDynamic ---

void test_bitvector_bitset_dynamic() {
  std::mt19937_64 rng(42);

  for (int64_t nbits : {64, 128, 200}) {
    auto vec = BitVector<BitsetDynamic>(10, nbits);
    REQUIRE(vec.size() == 10);
    REQUIRE(vec.nbits() == nbits);

    std::vector<BitsetDynamic> expected;
    for (int i = 0; i < 10; ++i) {
      BitsetDynamic val(nbits);
      for (int64_t b = 0; b < nbits; b += 64) {
        int64_t w = std::min(int64_t(64), nbits - b);
        val.set_range(b, w, rng()); // set_range masks to w bits internally
      }
      vec[i] = val;
      expected.push_back(val);
    }
    for (int i = 0; i < 10; ++i)
      REQUIRE(vec[i] == expected[i]);

    // const iteration
    int idx = 0;
    for (auto v : (BitVector<BitsetDynamic> const &)vec)
      REQUIRE(v == expected[idx++]);

    // operator== / !=
    auto v1 = BitVector<BitsetDynamic>(3, nbits);
    auto v2 = BitVector<BitsetDynamic>(3, nbits);
    for (int i = 0; i < 3; ++i) {
      v1[i] = expected[i];
      v2[i] = expected[i];
    }
    REQUIRE(v1 == v2);
    v2[1] = BitsetDynamic(nbits); // zero
    REQUIRE(v1 != v2);
  }
}

// --- BitsetStatic* (BitsetStatic1/2/4/8) ---

template <typename value_t> void test_bitvector_static_bitset() {
  constexpr int64_t nbits = std::numeric_limits<value_t>::digits;
  std::mt19937_64 rng(42);

  auto vec = BitVector<value_t>(10, nbits);
  REQUIRE(vec.size() == 10);
  REQUIRE(vec.nbits() == nbits);

  std::vector<value_t> expected;
  for (int i = 0; i < 10; ++i) {
    value_t val{};
    for (int64_t b = 0; b < nbits; b += 64)
      val.set_range(b, 64, rng());
    vec[i] = val;
    expected.push_back(val);
  }
  for (int i = 0; i < 10; ++i)
    REQUIRE(vec[i] == expected[i]);

  // const iteration
  int idx = 0;
  for (auto v : (BitVector<value_t> const &)vec)
    REQUIRE(v == expected[idx++]);

  // operator== / !=
  auto v1 = BitVector<value_t>(3, nbits);
  auto v2 = BitVector<value_t>(3, nbits);
  for (int i = 0; i < 3; ++i) {
    v1[i] = expected[i];
    v2[i] = expected[i];
  }
  REQUIRE(v1 == v2);
  v2[1] = value_t{};
  REQUIRE(v1 != v2);
}

// --- BitArray<uint64_t, N> (BitArray1..8) ---

template <typename value_t> void test_bitvector_bitarray() {
  constexpr int64_t nbits = 64; // raw uint64_t storage per element
  constexpr int N = value_t::nbits;
  constexpr int64_t slots = nbits / N;
  std::mt19937 rng(42);

  auto vec = BitVector<value_t>(10, nbits);
  REQUIRE(vec.size() == 10);
  REQUIRE(vec.nbits() == nbits);

  std::vector<std::vector<int>> exp(10);
  for (int i = 0; i < 10; ++i) {
    value_t elem = vec[i]; // zero element
    for (int64_t s = 0; s < slots; ++s) {
      int v =
          std::uniform_int_distribution<int>(0, (1 << N) - 1)(rng);
      elem.set(s, v);
      exp[i].push_back(v);
    }
    vec[i] = elem;
  }
  for (int i = 0; i < 10; ++i) {
    value_t elem = vec[i];
    for (int64_t s = 0; s < slots; ++s)
      REQUIRE(elem.get(s) == exp[i][s]);
  }

  // operator== / !=
  {
    auto v1 = BitVector<value_t>(3, nbits);
    auto v2 = BitVector<value_t>(3, nbits);
    REQUIRE(v1 == v2);
    value_t one{};
    one.set(0, 1);
    v1[0] = one;
    REQUIRE(v1 != v2);
  }
}

// --- BitArray<BitsetDynamic, N> (BitArrayLong1..8) ---

template <typename value_t> void test_bitvector_bitarray_long() {
  // Use 64-bit BitsetDynamic storage per element
  constexpr int64_t nbits = 64;
  constexpr int N = value_t::nbits;
  constexpr int64_t slots = nbits / N;
  std::mt19937 rng(42);

  auto vec = BitVector<value_t>(10, nbits);
  REQUIRE(vec.size() == 10);
  REQUIRE(vec.nbits() == nbits);

  std::vector<std::vector<int>> exp(10);
  for (int i = 0; i < 10; ++i) {
    value_t elem = vec[i]; // zero element with 64-bit BitsetDynamic storage
    for (int64_t s = 0; s < slots; ++s) {
      int v =
          std::uniform_int_distribution<int>(0, (1 << N) - 1)(rng);
      elem.set(s, v);
      exp[i].push_back(v);
    }
    vec[i] = elem;
  }
  for (int i = 0; i < 10; ++i) {
    value_t elem = vec[i];
    for (int64_t s = 0; s < slots; ++s)
      REQUIRE(elem.get(s) == exp[i][s]);
  }
}

TEST_CASE("bitvector", "[bits]") try {

  SECTION("integral") {
    Log("Testing bitvector - integral");
    test_bitvector_integral<uint16_t>();
    test_bitvector_integral<uint32_t>();
    test_bitvector_integral<uint64_t>();
  }

  SECTION("bitset_dynamic") {
    Log("Testing bitvector - bitset_dynamic");
    test_bitvector_bitset_dynamic();
  }

  SECTION("bitset_static") {
    Log("Testing bitvector - bitset_static");
    test_bitvector_static_bitset<BitsetStatic1>();
    test_bitvector_static_bitset<BitsetStatic2>();
    test_bitvector_static_bitset<BitsetStatic4>();
    test_bitvector_static_bitset<BitsetStatic8>();
  }

  SECTION("bitarray") {
    Log("Testing bitvector - bitarray");
    test_bitvector_bitarray<BitArray1>();
    test_bitvector_bitarray<BitArray2>();
    test_bitvector_bitarray<BitArray3>();
    test_bitvector_bitarray<BitArray4>();
    // Widths 5-7 are no longer compiled (boson dispatch promotes them to 8).
    test_bitvector_bitarray<BitArray8>();
  }

  SECTION("bitarray_long") {
    Log("Testing bitvector - bitarray_long");
    test_bitvector_bitarray_long<BitArrayLong1>();
    test_bitvector_bitarray_long<BitArrayLong2>();
    test_bitvector_bitarray_long<BitArrayLong3>();
    test_bitvector_bitarray_long<BitArrayLong4>();
    // Widths 5-7 are no longer compiled (boson dispatch promotes them to 8).
    test_bitvector_bitarray_long<BitArrayLong8>();
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
