// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../tests/catch.hpp"

#include <vector>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag::combinatorics;
using namespace xdiag::bits;
using namespace xdiag::math;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Decode the i-th digit of idx in base bound
static int64_t digit(int64_t idx, int64_t bound, int i) {
  for (int k = 0; k < i; ++k)
    idx /= bound;
  return idx % bound;
}

// ---------------------------------------------------------------------------
// Test: n(), bound(), size() accessors and iteration count
// ---------------------------------------------------------------------------
template <typename bitarray_t> void test_size(int64_t bound) {
  for (int n = 0; n <= 4; ++n) {
    BoundedMultisets<bitarray_t> ms(n, bound);
    REQUIRE(ms.n() == n);
    REQUIRE(ms.d() == bound);
    REQUIRE(ms.size() == ipow(bound, n));

    int64_t ctr = 0;
    for ([[maybe_unused]] auto seq : ms)
      ++ctr;
    REQUIRE(ctr == ms.size());
  }
}

// ---------------------------------------------------------------------------
// Test: every element of every sequence is in [0, bound)
// ---------------------------------------------------------------------------
template <typename bitarray_t> void test_element_bounds(int64_t bound) {
  for (int n = 1; n <= 4; ++n) {
    BoundedMultisets<bitarray_t> ms(n, bound);
    for (auto seq : ms)
      for (int i = 0; i < n; ++i) {
        REQUIRE(seq.get(i) >= 0);
        REQUIRE(seq.get(i) < bound);
      }
  }
}

// ---------------------------------------------------------------------------
// Test: exact ordering — sequence at linear index k encodes k in base bound
// (little-endian: slot 0 is the least significant digit)
// ---------------------------------------------------------------------------
template <typename bitarray_t> void test_ordering(int n, int64_t bound) {
  BoundedMultisets<bitarray_t> ms(n, bound);
  int64_t idx = 0;
  for (auto seq : ms) {
    for (int i = 0; i < n; ++i)
      REQUIRE(seq.get(i) == digit(idx, bound, i));
    ++idx;
  }
  REQUIRE(idx == ms.size());
}

// ---------------------------------------------------------------------------
// Test: n=1 — the single-element sequences are exactly 0, 1, …, bound-1
// ---------------------------------------------------------------------------
template <typename bitarray_t> void test_n1(int64_t bound) {
  BoundedMultisets<bitarray_t> ms(1, bound);
  REQUIRE(ms.size() == bound);
  int64_t expected = 0;
  for (auto seq : ms) {
    REQUIRE(seq.get(0) == expected);
    ++expected;
  }
}

// ---------------------------------------------------------------------------
// Test: n=0 — exactly one empty sequence
// ---------------------------------------------------------------------------
template <typename bitarray_t> void test_n0(int64_t bound) {
  BoundedMultisets<bitarray_t> ms(0, bound);
  REQUIRE(ms.size() == 1);
  int64_t ctr = 0;
  for ([[maybe_unused]] auto seq : ms)
    ++ctr;
  REQUIRE(ctr == 1);
}

// ---------------------------------------------------------------------------
// Test: all sequences are distinct
// ---------------------------------------------------------------------------
template <typename bitarray_t> void test_distinct(int n, int64_t bound) {
  BoundedMultisets<bitarray_t> ms(n, bound);
  std::vector<bitarray_t> seqs;
  seqs.reserve(ms.size());
  for (auto seq : ms)
    seqs.push_back(seq);

  // Adjacent pairs differ (ordering guarantees all-distinct)
  for (int64_t i = 1; i < (int64_t)seqs.size(); ++i)
    REQUIRE(seqs[i] != seqs[i - 1]);
}

// ---------------------------------------------------------------------------
// Test: the sequences are growing in size
// ---------------------------------------------------------------------------
template <typename bitarray_t> void test_growing(int n, int64_t bound) {
  BoundedMultisets<bitarray_t> ms(n, bound);
  bitarray_t prev;
  int64_t ctr = 0;
  for (auto seq : ms) {
    if (ctr > 0) {
      REQUIRE(seq > prev);
    }
    prev = seq;
    ++ctr;
  }
}

// ---------------------------------------------------------------------------
// Test: operator[] and index() — random access and round-trip
// ---------------------------------------------------------------------------
template <typename bitarray_t> void test_random_access(int n, int64_t bound) {
  BoundedMultisets<bitarray_t> ms(n, bound);

  // Collect elements sequentially
  std::vector<bitarray_t> elems;
  for (auto seq : ms)
    elems.push_back(seq);

  // operator[]: ms[i] matches sequential order
  for (int64_t i = 0; i < ms.size(); ++i)
    REQUIRE(ms[i] == elems[i]);

  // index: round-trip index(ms[i]) == i
  for (int64_t i = 0; i < ms.size(); ++i)
    REQUIRE(ms.index(elems[i]) == i);
}

// ---------------------------------------------------------------------------
// Test: iterator operator+ and operator+=
// ---------------------------------------------------------------------------
template <typename bitarray_t>
void test_iterator_advance(int n, int64_t bound) {
  BoundedMultisets<bitarray_t> ms(n, bound);
  if (ms.size() == 0)
    return;

  std::vector<bitarray_t> elems;
  for (auto seq : ms)
    elems.push_back(seq);

  // operator+: begin() + i matches elems[i]
  for (int64_t i = 0; i < ms.size(); ++i)
    REQUIRE(*(ms.begin() + i) == elems[i]);

  // operator+=: step-by-step advance
  auto it = ms.begin();
  for (int64_t i = 0; i < ms.size() - 1; ++i) {
    it += 1;
    REQUIRE(*it == elems[i + 1]);
  }

  // operator+=: larger jump
  if (ms.size() >= 4) {
    auto it2 = ms.begin();
    it2 += 3;
    REQUIRE(*it2 == elems[3]);
  }
}

// ---------------------------------------------------------------------------
// Test: equality operator
// ---------------------------------------------------------------------------
template <typename bitarray_t> void test_equality() {
  BoundedMultisets<bitarray_t> a(3, 2), b(3, 2), c(4, 2), d(3, 3);
  REQUIRE(a == b);
  REQUIRE(a != c); // different n
  REQUIRE(a != d); // different bound
}

// ---------------------------------------------------------------------------
// TEST_CASE
// ---------------------------------------------------------------------------

TEST_CASE("BoundedMultisets", "[combinatorics/bounded_multisets]") {

  SECTION("size and iteration count — uint64_t") {
    xdiag::Log(
        "Testing BoundedMultisets - size and iteration count — uint64_t");
    using A2 = BitArray<uint64_t, 1>;
    using A3 = BitArray<uint64_t, 2>;
    using A5 = BitArray<uint64_t, 3>;
    test_size<A2>(2);
    test_size<A3>(3);
    test_size<A3>(4);
    test_size<A5>(5);
    test_size<A5>(8);
  }

  SECTION("size and iteration count — uint32_t") {
    xdiag::Log(
        "Testing BoundedMultisets - size and iteration count — uint32_t");
    test_size<BitArray<uint32_t, 2>>(3);
  }

  SECTION("size and iteration count — Bitset types") {
    xdiag::Log(
        "Testing BoundedMultisets - size and iteration count — Bitset types");
    test_size<BitArray<BitsetStatic1, 1>>(2);
    test_size<BitArray<BitsetStatic2, 2>>(3);
    test_size<BitArray<BitsetStatic4, 3>>(5);
  }

  SECTION("element values in [0, bound)") {
    xdiag::Log("Testing BoundedMultisets - element values in [0, bound)");
    test_element_bounds<BitArray<uint64_t, 2>>(3);
    test_element_bounds<BitArray<uint64_t, 3>>(5);
    test_element_bounds<BitArray<uint32_t, 2>>(4);
  }

  SECTION("exact ordering") {
    xdiag::Log("Testing BoundedMultisets - exact ordering");
    test_ordering<BitArray<uint64_t, 2>>(3, 3);
    test_ordering<BitArray<uint64_t, 3>>(4, 5);
    test_ordering<BitArray<uint64_t, 1>>(8, 2);
    test_ordering<BitArray<uint32_t, 2>>(3, 4);
    test_ordering<BitArray<BitsetStatic1, 2>>(3, 3);
  }

  SECTION("n=1 yields 0..bound-1 in order") {
    xdiag::Log("Testing BoundedMultisets - n=1 yields 0..bound-1 in order");
    test_n1<BitArray<uint64_t, 2>>(3);
    test_n1<BitArray<uint64_t, 3>>(8);
  }

  SECTION("n=0 yields one empty sequence") {
    xdiag::Log("Testing BoundedMultisets - n=0 yields one empty sequence");
    test_n0<BitArray<uint64_t, 2>>(2);
    test_n0<BitArray<uint64_t, 3>>(5);
    test_n0<BitArray<uint32_t, 2>>(3);
  }

  SECTION("all sequences are distinct") {
    xdiag::Log("Testing BoundedMultisets - all sequences are distinct");
    test_distinct<BitArray<uint64_t, 2>>(3, 3);
    test_distinct<BitArray<uint64_t, 3>>(3, 5);
    test_distinct<BitArray<uint32_t, 2>>(4, 4);
  }

  SECTION("all sequences are growing in size") {
    xdiag::Log("Testing BoundedMultisets - all sequences are growing in size");
    test_growing<BitArray<uint64_t, 2>>(3, 3);
    test_growing<BitArray<uint64_t, 3>>(3, 5);
    test_growing<BitArray<uint32_t, 2>>(4, 4);
  }

  SECTION("equality operator") {
    xdiag::Log("Testing BoundedMultisets - equality operator");
    test_equality<BitArray<uint64_t, 2>>();
    test_equality<BitArray<uint32_t, 2>>();
    test_equality<BitArray<BitsetStatic2, 2>>();
  }

  SECTION("random access and index") {
    xdiag::Log("Testing BoundedMultisets - random access and index");
    test_random_access<BitArray<uint64_t, 2>>(3, 3);
    test_random_access<BitArray<uint64_t, 3>>(4, 5);
    test_random_access<BitArray<uint64_t, 1>>(6, 2);
    test_random_access<BitArray<uint32_t, 2>>(3, 4);
    test_random_access<BitArray<BitsetStatic2, 2>>(3, 3);
  }

  SECTION("iterator advance (+ and +=)") {
    xdiag::Log("Testing BoundedMultisets - iterator advance (+ and +=)");
    test_iterator_advance<BitArray<uint64_t, 2>>(3, 3);
    test_iterator_advance<BitArray<uint64_t, 3>>(4, 5);
    test_iterator_advance<BitArray<uint64_t, 1>>(6, 2);
    test_iterator_advance<BitArray<uint32_t, 2>>(3, 4);
  }
}
