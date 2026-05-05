// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../tests/catch.hpp"

#include <vector>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag::combinatorics;
using namespace xdiag::bits;
using namespace xdiag::math;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Brute-force reference: count all sequences from BoundedMultisets that sum to
// total, for cross-checking size()
template <int64_t nbits_val>
static int64_t brute_force_count(int n, int64_t total, int64_t bound) {
  using A = BitArray<uint64_t, nbits_val>;
  int64_t count = 0;
  int64_t total_seqs = ipow(bound, n);
  for (int64_t idx = 0; idx < total_seqs; ++idx) {
    int64_t tmp = idx, sum = 0;
    for (int i = 0; i < n; ++i) {
      sum += tmp % bound;
      tmp /= bound;
    }
    if (sum == total)
      ++count;
  }
  return count;
}

// ---------------------------------------------------------------------------
// Test: accessors and size matches brute-force count
// ---------------------------------------------------------------------------
template <typename bitarray_t, int64_t nbits_val>
void test_size(int n, int64_t total, int64_t bound) {
  BoundedPartitions<bitarray_t> bp(n, total, bound);
  REQUIRE(bp.n() == n);
  REQUIRE(bp.total() == total);
  REQUIRE(bp.bound() == bound);
  REQUIRE(bp.size() == brute_force_count<nbits_val>(n, total, bound));
}

// ---------------------------------------------------------------------------
// Test: iteration count matches size(), sum constraint, element bounds
// ---------------------------------------------------------------------------
template <typename bitarray_t>
void test_iter(int n, int64_t total, int64_t bound) {
  BoundedPartitions<bitarray_t> bp(n, total, bound);
  int64_t ctr = 0;
  for (auto seq : bp) {
    int64_t sum = 0;
    for (int i = 0; i < n; ++i) {
      REQUIRE(seq.get(i) >= 0);
      REQUIRE(seq.get(i) < bound);
      sum += seq.get(i);
    }
    REQUIRE(sum == total);
    ++ctr;
  }
  REQUIRE(ctr == bp.size());
}

// ---------------------------------------------------------------------------
// Test: sequences are in reverse-lexicographic order and all distinct
// ---------------------------------------------------------------------------
template <typename bitarray_t>
void test_order(int n, int64_t total, int64_t bound) {
  BoundedPartitions<bitarray_t> bp(n, total, bound);
  std::vector<bitarray_t> seqs;
  for (auto seq : bp)
    seqs.push_back(seq);

  for (int64_t k = 1; k < (int64_t)seqs.size(); ++k) {
    // Scan from most-significant slot (n-1) downward; first differing slot
    // must be strictly greater in seqs[k] (rlex order)
    bool rlex_greater = false;
    for (int i = n - 1; i >= 0; --i) {
      if (seqs[k].get(i) > seqs[k - 1].get(i)) {
        rlex_greater = true;
        break;
      }
      if (seqs[k].get(i) < seqs[k - 1].get(i))
        break;
    }
    REQUIRE(rlex_greater);
  }
}

// ---------------------------------------------------------------------------
// Test: exact sequence contents for n=3, bound=3, total=4
// Expected (rlex order): (2,2,0),(2,1,1),(1,2,1),(2,0,2),(1,1,2),(0,2,2)
// ---------------------------------------------------------------------------
void test_exact() {
  using A = BitArray<uint64_t, 2>;
  BoundedPartitions<A> bp(3, 4, 3);
  REQUIRE(bp.size() == 6);

  auto it = bp.begin();
  auto check = [&](int64_t a0, int64_t a1, int64_t a2) {
    REQUIRE((*it).get(0) == a0);
    REQUIRE((*it).get(1) == a1);
    REQUIRE((*it).get(2) == a2);
    ++it;
  };
  check(2, 2, 0);
  check(2, 1, 1);
  check(1, 2, 1);
  check(2, 0, 2);
  check(1, 1, 2);
  check(0, 2, 2);
  REQUIRE(it == bp.end());
}

// ---------------------------------------------------------------------------
// Test: n=0, total=0 -> one empty sequence; n=0, total>0 -> zero sequences
// ---------------------------------------------------------------------------
void test_n0() {
  using A = BitArray<uint64_t, 2>;
  BoundedPartitions<A> bp0(0, 0, 3);
  REQUIRE(bp0.size() == 1);
  int64_t ctr = 0;
  for ([[maybe_unused]] auto seq : bp0)
    ++ctr;
  REQUIRE(ctr == 1);

  BoundedPartitions<A> bp1(0, 1, 3);
  REQUIRE(bp1.size() == 0);
  ctr = 0;
  for ([[maybe_unused]] auto seq : bp1)
    ++ctr;
  REQUIRE(ctr == 0);
}

// ---------------------------------------------------------------------------
// Test: n=1 -> single-element sequences; only (total,) if total < bound
// ---------------------------------------------------------------------------
void test_n1() {
  using A = BitArray<uint64_t, 2>;
  BoundedPartitions<A> bp(1, 2, 3);
  REQUIRE(bp.size() == 1);
  auto it = bp.begin();
  REQUIRE((*it).get(0) == 2);
  REQUIRE(++it == bp.end());

  // total >= bound -> no sequences
  BoundedPartitions<A> bp2(1, 3, 3);
  REQUIRE(bp2.size() == 0);
}

// ---------------------------------------------------------------------------
// Test: total=0 -> only the all-zeros sequence
// ---------------------------------------------------------------------------
void test_total0() {
  using A = BitArray<uint64_t, 2>;
  for (int n = 1; n <= 5; ++n) {
    BoundedPartitions<A> bp(n, 0, 3);
    REQUIRE(bp.size() == 1);
    auto seq = *bp.begin();
    for (int i = 0; i < n; ++i)
      REQUIRE(seq.get(i) == 0);
  }
}

// ---------------------------------------------------------------------------
// Test: total = n*(bound-1) -> only the all-(bound-1) sequence
// ---------------------------------------------------------------------------
void test_total_max() {
  using A = BitArray<uint64_t, 2>;
  for (int n = 1; n <= 5; ++n) {
    int64_t bound = 3, total = n * (bound - 1);
    BoundedPartitions<A> bp(n, total, bound);
    REQUIRE(bp.size() == 1);
    auto seq = *bp.begin();
    for (int i = 0; i < n; ++i)
      REQUIRE(seq.get(i) == bound - 1);
  }
}

// ---------------------------------------------------------------------------
// Test: the sequences are growing in size
// ---------------------------------------------------------------------------
template <typename bitarray_t>
void test_growing(int n, int total, int64_t bound) {
  BoundedPartitions<bitarray_t> ms(n, total, bound);
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
template <typename bitarray_t>
void test_random_access(int n, int64_t total, int64_t bound) {
  BoundedPartitions<bitarray_t> bp(n, total, bound);

  std::vector<bitarray_t> elems;
  for (auto seq : bp)
    elems.push_back(seq);

  // operator[]: bp[i] matches sequential order
  for (int64_t i = 0; i < bp.size(); ++i)
    REQUIRE(bp[i] == elems[i]);

  // index: round-trip index(bp[i]) == i
  for (int64_t i = 0; i < bp.size(); ++i)
    REQUIRE(bp.index(elems[i]) == i);
}

// ---------------------------------------------------------------------------
// Test: iterator operator+ and operator+=
// ---------------------------------------------------------------------------
template <typename bitarray_t>
void test_iterator_advance(int n, int64_t total, int64_t bound) {
  BoundedPartitions<bitarray_t> bp(n, total, bound);
  if (bp.size() == 0)
    return;

  std::vector<bitarray_t> elems;
  for (auto seq : bp)
    elems.push_back(seq);

  // operator+: begin() + i matches elems[i]
  for (int64_t i = 0; i < bp.size(); ++i)
    REQUIRE(*(bp.begin() + i) == elems[i]);

  // operator+=: step-by-step advance
  auto it = bp.begin();
  for (int64_t i = 0; i < bp.size() - 1; ++i) {
    it += 1;
    REQUIRE(*it == elems[i + 1]);
  }

  // operator+=: larger jump
  if (bp.size() >= 4) {
    auto it2 = bp.begin();
    it2 += 3;
    REQUIRE(*it2 == elems[3]);
  }
}

// ---------------------------------------------------------------------------
// Test: equality operator
// ---------------------------------------------------------------------------
void test_equality() {
  using A = BitArray<uint64_t, 2>;
  BoundedPartitions<A> a(3, 4, 3), b(3, 4, 3), c(3, 3, 3), d(3, 4, 4),
      e(4, 4, 3);
  REQUIRE(a == b);
  REQUIRE(a != c); // different total
  REQUIRE(a != d); // different bound
  REQUIRE(a != e); // different n
}

// ---------------------------------------------------------------------------
// TEST_CASE
// ---------------------------------------------------------------------------

TEST_CASE("BoundedPartitions", "[combinatorics/bounded_partitions]") {

  SECTION("size matches brute-force count") {
    xdiag::Log("Testing BoundedPartitions - size matches brute-force count");
    // Use template param nbits_val=2 for brute_force_count (bound <= 4)
    using A2 = BitArray<uint64_t, 2>;
    for (int n = 0; n <= 4; ++n)
      for (int64_t t = 0; t <= n * 2; ++t)
        test_size<A2, 2>(n, t, 3);

    // bound=2 with nbits=1
    using A1 = BitArray<uint64_t, 1>;
    for (int n = 0; n <= 6; ++n)
      for (int64_t t = 0; t <= n; ++t)
        test_size<A1, 1>(n, t, 2);
  }

  SECTION("iteration: count, sum constraint, element bounds") {
    xdiag::Log("Testing BoundedPartitions - iteration: count, sum constraint, "
               "element bounds");
    test_iter<BitArray<uint64_t, 2>>(3, 4, 3);
    test_iter<BitArray<uint64_t, 2>>(4, 5, 4);
    test_iter<BitArray<uint64_t, 3>>(5, 7, 5);
    test_iter<BitArray<uint32_t, 2>>(3, 3, 3);
    test_iter<BitArray<BitsetStatic1, 2>>(3, 4, 3);
  }

  SECTION("reverse-lexicographic order") {
    xdiag::Log("Testing BoundedPartitions - reverse-lexicographic order");
    test_order<BitArray<uint64_t, 2>>(3, 4, 3);
    test_order<BitArray<uint64_t, 3>>(4, 6, 5);
    test_order<BitArray<uint64_t, 1>>(6, 3, 2);
    test_order<BitArray<uint32_t, 2>>(4, 5, 4);
  }

  SECTION("exact sequence order (n=3, bound=3, total=4)") {
    xdiag::Log("Testing BoundedPartitions - exact sequence order (n=3, "
               "bound=3, total=4)");
    test_exact();
  }

  SECTION("n=0 edge cases") {
    xdiag::Log("Testing BoundedPartitions - n=0 edge cases");
    test_n0();
  }

  SECTION("n=1 edge cases") {
    xdiag::Log("Testing BoundedPartitions - n=1 edge cases");
    test_n1();
  }

  SECTION("total=0 yields all-zeros") {
    xdiag::Log("Testing BoundedPartitions - total=0 yields all-zeros");
    test_total0();
  }

  SECTION("total=n*(bound-1) yields all-(bound-1)") {
    xdiag::Log(
        "Testing BoundedPartitions - total=n*(bound-1) yields all-(bound-1)");
    test_total_max();
  }

  SECTION("all sequences are growing in size") {
    xdiag::Log("Testing BoundedPartitions - all sequences are growing in size");
    test_growing<BitArray<uint64_t, 2>>(3, 3, 3);
    test_growing<BitArray<uint64_t, 3>>(3, 5, 5);
    test_growing<BitArray<uint32_t, 2>>(4, 4, 4);
  }

  SECTION("equality operator") {
    xdiag::Log("Testing BoundedPartitions - equality operator");
    test_equality();
  }

  SECTION("random access and index") {
    xdiag::Log("Testing BoundedPartitions - random access and index");
    test_random_access<BitArray<uint64_t, 2>>(3, 4, 3);
    test_random_access<BitArray<uint64_t, 2>>(4, 5, 4);
    test_random_access<BitArray<uint64_t, 3>>(5, 7, 5);
    test_random_access<BitArray<uint64_t, 1>>(6, 3, 2);
    test_random_access<BitArray<uint32_t, 2>>(3, 3, 3);
    test_random_access<BitArray<BitsetStatic1, 2>>(3, 4, 3);
  }

  SECTION("iterator advance (+ and +=)") {
    xdiag::Log("Testing BoundedPartitions - iterator advance (+ and +=)");
    test_iterator_advance<BitArray<uint64_t, 2>>(3, 4, 3);
    test_iterator_advance<BitArray<uint64_t, 2>>(4, 5, 4);
    test_iterator_advance<BitArray<uint64_t, 3>>(5, 7, 5);
    test_iterator_advance<BitArray<uint32_t, 2>>(3, 3, 3);
  }
}
