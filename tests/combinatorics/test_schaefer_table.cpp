// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <vector>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/utils/logger.hpp>

// Iterate SchaeferTable and verify count, per-element sum, and value bounds.
template <typename bitarray_t>
void test_schaefer_iteration(int64_t bound) {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 0; n <= 6; ++n) {
    for (int total = 0; total <= n * (bound - 1); ++total) {
      SchaeferTable<bitarray_t> st(n, total, bound);
      REQUIRE(st.n() == n);
      REQUIRE(st.total() == total);
      REQUIRE(st.d() == bound);

      int64_t ctr = 0;
      for (auto seq : st) {
        int64_t sum = 0;
        for (int i = 0; i < n; ++i) {
          REQUIRE(seq.get(i) >= 0);
          REQUIRE(seq.get(i) < bound);
          sum += seq.get(i);
        }
        REQUIRE(sum == total);
        ++ctr;
      }
      REQUIRE(ctr == st.size());
    }
  }
}

// index() must equal the sequential position in iteration.
template <typename bitarray_t>
void test_schaefer_index(int64_t bound) {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 0; n <= 6; ++n) {
    for (int total = 0; total <= n * (bound - 1); ++total) {
      SchaeferTable<bitarray_t> st(n, total, bound);
      int64_t idx = 0;
      for (auto seq : st) {
        REQUIRE(st.index(seq) == idx);
        ++idx;
      }
    }
  }
}

// operator[](i) matches sequential order; index() is its inverse.
template <typename bitarray_t>
void test_schaefer_random_access(int64_t bound) {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 0; n <= 6; ++n) {
    for (int total = 0; total <= n * (bound - 1); ++total) {
      SchaeferTable<bitarray_t> st(n, total, bound);

      std::vector<bitarray_t> elems;
      for (auto seq : st)
        elems.push_back(seq);

      for (int64_t i = 0; i < st.size(); ++i)
        REQUIRE(st[i] == elems[i]);

      for (int64_t i = 0; i < st.size(); ++i)
        REQUIRE(st.index(elems[i]) == i);
    }
  }
}

// Iterator operator+ and operator+=
template <typename bitarray_t>
void test_schaefer_iterator_advance(int64_t bound) {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 1; n <= 6; ++n) {
    for (int total = 0; total <= n * (bound - 1); ++total) {
      SchaeferTable<bitarray_t> st(n, total, bound);
      if (st.size() == 0)
        continue;

      std::vector<bitarray_t> elems;
      for (auto seq : st)
        elems.push_back(seq);

      for (int64_t i = 0; i < st.size(); ++i)
        REQUIRE(*(st.begin() + i) == elems[i]);

      auto it = st.begin();
      for (int64_t i = 0; i < st.size() - 1; ++i) {
        it += 1;
        REQUIRE(*it == elems[i + 1]);
      }

      if (st.size() >= 3) {
        auto it2 = st.begin();
        it2 += 2;
        REQUIRE(*it2 == elems[2]);
      }
    }
  }
}

// SchaeferTable must produce the exact same sequence as BoundedPartitions.
template <typename bitarray_t>
void test_schaefer_vs_bounded_partitions(int64_t bound) {
  using namespace xdiag;
  using namespace xdiag::combinatorics;

  for (int n = 0; n <= 6; ++n) {
    for (int total = 0; total <= n * (bound - 1); ++total) {
      SchaeferTable<bitarray_t> st(n, total, bound);
      BoundedPartitions<bitarray_t> bp(n, total, bound);

      REQUIRE(st.n() == bp.n());
      REQUIRE(st.total() == bp.total());
      REQUIRE(st.d() == bp.d());
      REQUIRE(st.size() == bp.size());

      auto it_st = st.begin();
      auto it_bp = bp.begin();
      while (it_st != st.end()) {
        REQUIRE(*it_st == *it_bp);
        ++it_st;
        ++it_bp;
      }
      REQUIRE(it_bp == bp.end());
    }
  }
}

TEST_CASE("SchaeferTable", "[combinatorics]") {
  using namespace xdiag::bits;
  using A2u64 = BitArray<uint64_t, 2>; // 2-bit slots: bound up to 4
  using A3u64 = BitArray<uint64_t, 3>; // 3-bit slots: bound up to 8
  using A2u32 = BitArray<uint32_t, 2>; // smaller storage


  SECTION("iteration (bound=2)") {
    xdiag::Log("Testing SchaeferTable - iteration (bound=2)");
    test_schaefer_iteration<A2u64>(2);
    test_schaefer_iteration<A3u64>(2);
    test_schaefer_iteration<A2u32>(2);
  }

  SECTION("iteration (bound=3)") {
    xdiag::Log("Testing SchaeferTable - iteration (bound=3)");
    test_schaefer_iteration<A2u64>(3);
    test_schaefer_iteration<A3u64>(3);
  }

  SECTION("iteration (bound=4)") {
    xdiag::Log("Testing SchaeferTable - iteration (bound=4)");
    test_schaefer_iteration<A2u64>(4);
    test_schaefer_iteration<A3u64>(4);
  }

  SECTION("index (bound=2)") {
    xdiag::Log("Testing SchaeferTable - index (bound=2)");
    test_schaefer_index<A2u64>(2);
    test_schaefer_index<A3u64>(2);
    test_schaefer_index<A2u32>(2);
  }

  SECTION("index (bound=3)") {
    xdiag::Log("Testing SchaeferTable - index (bound=3)");
    test_schaefer_index<A2u64>(3);
    test_schaefer_index<A3u64>(3);
  }

  SECTION("index (bound=4)") {
    xdiag::Log("Testing SchaeferTable - index (bound=4)");
    test_schaefer_index<A2u64>(4);
    test_schaefer_index<A3u64>(4);
  }

  SECTION("random access and index round-trip (bound=2)") {
    xdiag::Log("Testing SchaeferTable - random access and index round-trip (bound=2)");
    test_schaefer_random_access<A2u64>(2);
    test_schaefer_random_access<A3u64>(2);
  }

  SECTION("random access and index round-trip (bound=3)") {
    xdiag::Log("Testing SchaeferTable - random access and index round-trip (bound=3)");
    test_schaefer_random_access<A2u64>(3);
    test_schaefer_random_access<A3u64>(3);
  }

  SECTION("iterator advance (+ and +=) (bound=2)") {
    xdiag::Log("Testing SchaeferTable - iterator advance (+ and +=) (bound=2)");
    test_schaefer_iterator_advance<A2u64>(2);
    test_schaefer_iterator_advance<A3u64>(2);
  }

  SECTION("iterator advance (+ and +=) (bound=3)") {
    xdiag::Log("Testing SchaeferTable - iterator advance (+ and +=) (bound=3)");
    test_schaefer_iterator_advance<A2u64>(3);
    test_schaefer_iterator_advance<A3u64>(3);
  }

  SECTION("vs BoundedPartitions (bound=2)") {
    xdiag::Log("Testing SchaeferTable - vs BoundedPartitions (bound=2)");
    test_schaefer_vs_bounded_partitions<A2u64>(2);
    test_schaefer_vs_bounded_partitions<A3u64>(2);
    test_schaefer_vs_bounded_partitions<A2u32>(2);
  }

  SECTION("vs BoundedPartitions (bound=3)") {
    xdiag::Log("Testing SchaeferTable - vs BoundedPartitions (bound=3)");
    test_schaefer_vs_bounded_partitions<A2u64>(3);
    test_schaefer_vs_bounded_partitions<A3u64>(3);
  }

  SECTION("vs BoundedPartitions (bound=4)") {
    xdiag::Log("Testing SchaeferTable - vs BoundedPartitions (bound=4)");
    test_schaefer_vs_bounded_partitions<A2u64>(4);
    test_schaefer_vs_bounded_partitions<A3u64>(4);
  }
}
