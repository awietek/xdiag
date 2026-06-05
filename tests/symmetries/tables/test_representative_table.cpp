// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <cmath>
#include <vector>

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/config.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/symmetries/action/isrepresentative.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/tables/representative_table.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::bits;
using namespace xdiag::combinatorics;
using namespace xdiag::symmetries;

// -----------------------------------------------------------------------
// Generic invariant checker. The SitePermutation is derived from
// irrep.group() — identical to what the constructor uses internally.
//
// Checks:
//   1. size() in [0, enumeration.size()]
//   2. Iteration yields exactly size() elements, all isrepresentative
//   3. operator[](ri) matches iteration order
//   4. For every rep at index ri, every group element sym:
//        state = sp.apply(sym, rep)
//        representative(idx) == rep
//        representative_index(idx) == ri              (0-based)
//        raw_representative_index(idx) == ri + 1      (1-based storage)
//        sp.apply(representative_symmetry(idx), state) == rep
//   5. All norms are positive
//   6. table[representative_index(idx)] == representative(idx)
//      for every enumeration state with nonzero norm
//   7. Two identically-constructed tables compare equal
// -----------------------------------------------------------------------
template <typename enumeration_t>
void test_invariants(enumeration_t const &enumeration,
                     Representation const &irrep) {
  using bit_t = typename enumeration_t::bit_t;
  auto sp = SitePermutation(irrep.group());
  RepresentativeTable<enumeration_t> table(enumeration, irrep.group(), irrep.characters());

  REQUIRE(table.size() >= 0);
  REQUIRE(table.size() <= enumeration.size());

  // Iteration count and isrepresentative
  int64_t count = 0;
  for (auto rep : table) {
    REQUIRE(isrepresentative(rep, sp));
    ++count;
  }
  REQUIRE(count == table.size());

  // operator[] matches iteration order
  {
    int64_t ri = 0;
    for (auto rep : table) {
      REQUIRE(table[ri] == rep);
      ++ri;
    }
  }

  // Orbit structure for every representative
  for (int64_t ri = 0; ri < table.size(); ++ri) {
    bit_t rep = table[ri];
    for (int64_t sym = 0; sym < sp.size(); ++sym) {
      bit_t state = sp.apply(sym, rep);
      int64_t idx = enumeration.index(state);
      REQUIRE(table.representative(idx) == rep);
      REQUIRE(table.representative_index(idx) == ri);
      REQUIRE(table.raw_representative_index(idx) == ri + 1);
      REQUIRE(sp.apply(table.representative_symmetry(idx), state) == rep);
    }
  }

  // Norms are positive for all representatives
  for (int64_t ri = 0; ri < table.size(); ++ri) {
    REQUIRE(table.representative_norm(ri) > 0.0);
  }

  // Consistency: table[representative_index(idx)] == representative(idx)
  // for every enumeration state that has a nonzero norm (raw_index != 0).
  for (auto state : enumeration) {
    int64_t idx = enumeration.index(state);
    if (table.raw_representative_index(idx) != 0) {
      int64_t ri = table.representative_index(idx);
      REQUIRE(ri >= 0);
      REQUIRE(ri < table.size());
      REQUIRE(table[ri] == table.representative(idx));
    }
  }

  // Two tables from the same inputs compare equal
  RepresentativeTable<enumeration_t> table2(enumeration, irrep.group(), irrep.characters());
  REQUIRE(table == table2);
  REQUIRE_FALSE(table != table2);
}

template <typename enumeration_t>
void check_size(enumeration_t const &enumeration, Representation const &irrep,
                int64_t expected) {
  RepresentativeTable<enumeration_t> table(enumeration, irrep.group(), irrep.characters());
  REQUIRE(table.size() == expected);
}

// -----------------------------------------------------------------------
// TEST CASE
// -----------------------------------------------------------------------

TEST_CASE("representative_table", "[symmetries]") try {

  // Groups and irreps used throughout.
  // cyclic_group(4): rotations r0=id, r1, r2, r3 with chi_k(rm) = exp(2πikm/4)
  //   k=0: chi = {1, 1, 1, 1}         (trivial, real)
  //   k=1: chi = {1, i, -1, -i}       (complex)
  //   k=2: chi = {1, -1, 1, -1}       (real, non-trivial)
  auto irrep4_0 = cyclic_group_irrep(4, 0);
  auto irrep4_1 = cyclic_group_irrep(4, 1);
  auto irrep4_2 = cyclic_group_irrep(4, 2);

  auto irrep3_0 = cyclic_group_irrep(3, 0);
  auto irrep3_1 = cyclic_group_irrep(3, 1);

  // =====================================================================
  // Subsets — all 16 bit-strings of length 4, under cyclic_group(4)
  //
  // Orbits:
  //   {0}           size 1, stab=G  → norm=sqrt(4)=2
  //   {1,2,4,8}     size 4, stab={id} → norm=1
  //   {3,6,12,9}    size 4, stab={id} → norm=1
  //   {5,10}        size 2, stab={id,r2} → norm=sqrt(2)
  //   {7,14,13,11}  size 4, stab={id} → norm=1
  //   {15}          size 1, stab=G  → norm=sqrt(4)=2
  //
  // k=0: 6 reps (all orbits have nonzero norm)
  // k=2: 4 reps ({0} and {15} excluded: sum_stab chi_2 = 1+(-1)+1+(-1) = 0)
  // k=1: 3 reps ({0}, {5,10}, {15} excluded)
  // =====================================================================

  SECTION("Subsets<uint32_t>") {
    Log("Testing RepresentativeTable - Subsets<uint32_t>");
    auto subsets = Subsets<uint32_t>(4);
    auto sp4 = SitePermutation(irrep4_0.group());

    test_invariants(subsets, irrep4_0);
    test_invariants(subsets, irrep4_2);
    test_invariants(subsets, irrep4_1);
    check_size(subsets, irrep4_0, 6);
    check_size(subsets, irrep4_2, 4);
    check_size(subsets, irrep4_1, 3);

    // ---- Specific values: k=0 ----
    {
      RepresentativeTable<Subsets<uint32_t>> table(subsets, irrep4_0.group(), irrep4_0.characters());

      // Representatives in enumeration (lexicographic) order
      std::vector<uint32_t> expected_reps = {0, 1, 3, 5, 7, 15};
      std::vector<uint32_t> got_reps;
      for (auto r : table)
        got_reps.push_back(r);
      REQUIRE(got_reps == expected_reps);

      // operator[] matches
      for (int64_t ri = 0; ri < (int64_t)expected_reps.size(); ++ri)
        REQUIRE(table[ri] == expected_reps[ri]);

      // Norms by representative index
      REQUIRE(table.representative_norm(0) == Approx(2.0));
      REQUIRE(table.representative_norm(1) == Approx(1.0));
      REQUIRE(table.representative_norm(2) == Approx(1.0));
      REQUIRE(table.representative_norm(3) == Approx(std::sqrt(2.0)));
      REQUIRE(table.representative_norm(4) == Approx(1.0));
      REQUIRE(table.representative_norm(5) == Approx(2.0));

      // Orbit {1,2,4,8}: all map to rep 1 (ri=1), raw_index = 2
      for (uint32_t s : {1u, 2u, 4u, 8u}) {
        int64_t idx = subsets.index(s);
        REQUIRE(table.representative(idx) == 1u);
        REQUIRE(table.representative_index(idx) == 1);
        REQUIRE(table.raw_representative_index(idx) == 2);
        REQUIRE(sp4.apply(table.representative_symmetry(idx), s) == 1u);
      }

      // Orbit {5,10}: rep=5 (ri=3), raw_index = 4
      for (uint32_t s : {5u, 10u}) {
        int64_t idx = subsets.index(s);
        REQUIRE(table.representative(idx) == 5u);
        REQUIRE(table.representative_index(idx) == 3);
        REQUIRE(table.raw_representative_index(idx) == 4);
        REQUIRE(sp4.apply(table.representative_symmetry(idx), s) == 5u);
      }

      // Inequality: different irrep → tables differ
      RepresentativeTable<Subsets<uint32_t>> table_k2(subsets, irrep4_2.group(), irrep4_2.characters());
      REQUIRE(table != table_k2);
      REQUIRE_FALSE(table == table_k2);
    }

    // ---- Zero-norm state detection: k=2 ----
    // States 0 and 15 are fixed by all rotations; chi-sum = 0 → excluded.
    {
      RepresentativeTable<Subsets<uint32_t>> table(subsets, irrep4_2.group(), irrep4_2.characters());
      for (uint32_t s : {0u, 15u}) {
        int64_t idx = subsets.index(s);
        REQUIRE(table.raw_representative_index(idx) == 0);
        REQUIRE(table.representative_index(idx) == -1);
      }
    }
  }

  SECTION("Subsets<uint64_t>") {
    Log("Testing RepresentativeTable - Subsets<uint64_t>");
    auto subsets = Subsets<uint64_t>(4);
    test_invariants(subsets, irrep4_0);
    check_size(subsets, irrep4_0, 6);
  }

  // =====================================================================
  // Combinations(4,2) — 6 states under cyclic_group(4)
  //
  // Orbits: {3,6,12,9} (size 4) and {5,10} (size 2, stab={id,r2})
  //
  // k=0: 2 reps
  // k=2: 2 reps (chi_2(r2)=1, norm²({5,10})=1+1=2 > 0)
  // k=1: 1 rep  (chi_1(r2)=-1, norm²({5,10})=1+(-1)=0 → excluded)
  // =====================================================================

  SECTION("Combinations<uint32_t>") {
    Log("Testing RepresentativeTable - Combinations<uint32_t>");
    auto combos = Combinations<uint32_t>(4, 2);
    auto sp4 = SitePermutation(irrep4_0.group());

    test_invariants(combos, irrep4_0);
    test_invariants(combos, irrep4_2);
    test_invariants(combos, irrep4_1);
    check_size(combos, irrep4_0, 2);
    check_size(combos, irrep4_2, 2);
    check_size(combos, irrep4_1, 1);

    // Specific values for k=0
    {
      RepresentativeTable<Combinations<uint32_t>> table(combos, irrep4_0.group(), irrep4_0.characters());

      auto it = table.begin();
      REQUIRE(*it == 3u);
      ++it;
      REQUIRE(*it == 5u);
      ++it;
      REQUIRE(it == table.end());

      for (uint32_t s : {3u, 6u, 9u, 12u}) {
        int64_t idx = combos.index(s);
        REQUIRE(table.representative(idx) == 3u);
        REQUIRE(table.representative_index(idx) == 0);
        REQUIRE(sp4.apply(table.representative_symmetry(idx), s) == 3u);
      }
      for (uint32_t s : {5u, 10u}) {
        int64_t idx = combos.index(s);
        REQUIRE(table.representative(idx) == 5u);
        REQUIRE(table.representative_index(idx) == 1);
        REQUIRE(sp4.apply(table.representative_symmetry(idx), s) == 5u);
      }
    }

    // Zero-norm for k=1: orbit {5,10} excluded
    {
      RepresentativeTable<Combinations<uint32_t>> table(combos, irrep4_1.group(), irrep4_1.characters());
      for (uint32_t s : {5u, 10u}) {
        int64_t idx = combos.index(s);
        REQUIRE(table.raw_representative_index(idx) == 0);
        REQUIRE(table.representative_index(idx) == -1);
      }
    }
  }

  SECTION("Combinations<uint64_t>") {
    Log("Testing RepresentativeTable - Combinations<uint64_t>");
    auto combos = Combinations<uint64_t>(4, 2);
    test_invariants(combos, irrep4_0);
    check_size(combos, irrep4_0, 2);
  }

  // =====================================================================
  // Combinations (Bitset types) — same 6 states, same orbit structure
  // =====================================================================

  SECTION("Combinations<BitsetStatic2>") {
    Log("Testing RepresentativeTable - Combinations<BitsetStatic2>");
    auto combos = Combinations<BitsetStatic2>(4, 2);
    test_invariants(combos, irrep4_0);
    check_size(combos, irrep4_0, 2);
  }

  SECTION("Combinations<BitsetStatic4>") {
    Log("Testing RepresentativeTable - Combinations<BitsetStatic4>");
    auto combos = Combinations<BitsetStatic4>(4, 2);
    test_invariants(combos, irrep4_0);
    check_size(combos, irrep4_0, 2);
  }

  SECTION("Combinations<BitsetStatic8>") {
    Log("Testing RepresentativeTable - Combinations<BitsetStatic8>");
    auto combos = Combinations<BitsetStatic8>(4, 2);
    test_invariants(combos, irrep4_0);
    check_size(combos, irrep4_0, 2);
  }

  SECTION("Combinations<BitsetDynamic>") {
    Log("Testing RepresentativeTable - Combinations<BitsetDynamic>");
    auto combos = Combinations<BitsetDynamic>(4, 2);
    test_invariants(combos, irrep4_0);
    check_size(combos, irrep4_0, 2);
  }

  // =====================================================================
  // LinTable — O(1)-index drop-in for Combinations; same orbit structure
  // =====================================================================

  SECTION("LinTable<uint32_t>") {
    Log("Testing RepresentativeTable - LinTable<uint32_t>");
    auto lt = LinTable<uint32_t>(4, 2);
    test_invariants(lt, irrep4_0);
    test_invariants(lt, irrep4_2);
    check_size(lt, irrep4_0, 2);
  }

  SECTION("LinTable<uint64_t>") {
    Log("Testing RepresentativeTable - LinTable<uint64_t>");
    auto lt = LinTable<uint64_t>(4, 2);
    test_invariants(lt, irrep4_0);
    check_size(lt, irrep4_0, 2);
  }

  // =====================================================================
  // BoundedMultisets — with bound=2 equivalent to Subsets; same counts
  // =====================================================================

  SECTION("BoundedMultisets<BitArray1>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArray1>");
    auto ms = BoundedMultisets<BitArray1>(4, 2);
    test_invariants(ms, irrep4_0);
    test_invariants(ms, irrep4_2);
    check_size(ms, irrep4_0, 6);
    check_size(ms, irrep4_2, 4);
  }

  SECTION("BoundedMultisets<BitArray2>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArray2>");
    auto ms = BoundedMultisets<BitArray2>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray3>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArray3>");
    // bound=4, 3 bits per slot: 4^4 = 256 states
    auto ms = BoundedMultisets<BitArray3>(4, 4);
    test_invariants(ms, irrep4_0);
    test_invariants(ms, irrep4_1);
  }

  SECTION("BoundedMultisets<BitArray4>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArray4>");
    auto ms = BoundedMultisets<BitArray4>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray5>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArray5>");
    auto ms = BoundedMultisets<BitArray5>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray6>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArray6>");
    auto ms = BoundedMultisets<BitArray6>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray7>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArray7>");
    auto ms = BoundedMultisets<BitArray7>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray8>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArray8>");
    auto ms = BoundedMultisets<BitArray8>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong1>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArrayLong1>");
    auto ms = BoundedMultisets<BitArrayLong1>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
    check_size(ms, irrep4_2, 4);
  }

  SECTION("BoundedMultisets<BitArrayLong2>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArrayLong2>");
    auto ms = BoundedMultisets<BitArrayLong2>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong3>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArrayLong3>");
    auto ms = BoundedMultisets<BitArrayLong3>(4, 4);
    test_invariants(ms, irrep4_0);
    test_invariants(ms, irrep4_1);
  }

  SECTION("BoundedMultisets<BitArrayLong4>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArrayLong4>");
    auto ms = BoundedMultisets<BitArrayLong4>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong5>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArrayLong5>");
    auto ms = BoundedMultisets<BitArrayLong5>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong6>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArrayLong6>");
    auto ms = BoundedMultisets<BitArrayLong6>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong7>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArrayLong7>");
    auto ms = BoundedMultisets<BitArrayLong7>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong8>") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArrayLong8>");
    auto ms = BoundedMultisets<BitArrayLong8>(4, 2);
    test_invariants(ms, irrep4_0);
    check_size(ms, irrep4_0, 6);
  }

  // =====================================================================
  // BoundedPartitions
  //
  // n=4, total=2, bound=2 (BitArray1): same as Combinations(4,2) → 2 reps
  // n=4, total=2, bound=3 (BitArray2+): 10 states, 3 orbits → 3 reps
  // =====================================================================

  SECTION("BoundedPartitions<BitArray1>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArray1>");
    auto bp = BoundedPartitions<BitArray1>(4, 2, 2);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 2);
  }

  SECTION("BoundedPartitions<BitArray2>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArray2>");
    auto bp = BoundedPartitions<BitArray2>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    test_invariants(bp, irrep4_2);
    check_size(bp, irrep4_0, 3);
    check_size(bp, irrep4_2, 3);
  }

  SECTION("BoundedPartitions<BitArray3>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArray3>");
    auto bp = BoundedPartitions<BitArray3>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray4>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArray4>");
    auto bp = BoundedPartitions<BitArray4>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray5>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArray5>");
    auto bp = BoundedPartitions<BitArray5>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray6>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArray6>");
    auto bp = BoundedPartitions<BitArray6>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray7>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArray7>");
    auto bp = BoundedPartitions<BitArray7>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray8>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArray8>");
    auto bp = BoundedPartitions<BitArray8>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong1>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArrayLong1>");
    auto bp = BoundedPartitions<BitArrayLong1>(4, 2, 2);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 2);
  }

  SECTION("BoundedPartitions<BitArrayLong2>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArrayLong2>");
    auto bp = BoundedPartitions<BitArrayLong2>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong3>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArrayLong3>");
    auto bp = BoundedPartitions<BitArrayLong3>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong4>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArrayLong4>");
    auto bp = BoundedPartitions<BitArrayLong4>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong5>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArrayLong5>");
    auto bp = BoundedPartitions<BitArrayLong5>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong6>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArrayLong6>");
    auto bp = BoundedPartitions<BitArrayLong6>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong7>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArrayLong7>");
    auto bp = BoundedPartitions<BitArrayLong7>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong8>") {
    Log("Testing RepresentativeTable - BoundedPartitions<BitArrayLong8>");
    auto bp = BoundedPartitions<BitArrayLong8>(4, 2, 3);
    test_invariants(bp, irrep4_0);
    check_size(bp, irrep4_0, 3);
  }

  // =====================================================================
  // SchaeferTable — same enumeration as BoundedPartitions, O(1) index()
  // Only BitArray1..8
  // =====================================================================

  SECTION("SchaeferTable<BitArray1>") {
    Log("Testing RepresentativeTable - SchaeferTable<BitArray1>");
    auto st = SchaeferTable<BitArray1>(4, 2, 2);
    test_invariants(st, irrep4_0);
    check_size(st, irrep4_0, 2);
  }

  SECTION("SchaeferTable<BitArray2>") {
    Log("Testing RepresentativeTable - SchaeferTable<BitArray2>");
    auto st = SchaeferTable<BitArray2>(4, 2, 3);
    test_invariants(st, irrep4_0);
    test_invariants(st, irrep4_2);
    check_size(st, irrep4_0, 3);
    check_size(st, irrep4_2, 3);
  }

  SECTION("SchaeferTable<BitArray3>") {
    Log("Testing RepresentativeTable - SchaeferTable<BitArray3>");
    auto st = SchaeferTable<BitArray3>(4, 2, 3);
    test_invariants(st, irrep4_0);
    check_size(st, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray4>") {
    Log("Testing RepresentativeTable - SchaeferTable<BitArray4>");
    auto st = SchaeferTable<BitArray4>(4, 2, 3);
    test_invariants(st, irrep4_0);
    check_size(st, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray5>") {
    Log("Testing RepresentativeTable - SchaeferTable<BitArray5>");
    auto st = SchaeferTable<BitArray5>(4, 2, 3);
    test_invariants(st, irrep4_0);
    check_size(st, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray6>") {
    Log("Testing RepresentativeTable - SchaeferTable<BitArray6>");
    auto st = SchaeferTable<BitArray6>(4, 2, 3);
    test_invariants(st, irrep4_0);
    check_size(st, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray7>") {
    Log("Testing RepresentativeTable - SchaeferTable<BitArray7>");
    auto st = SchaeferTable<BitArray7>(4, 2, 3);
    test_invariants(st, irrep4_0);
    check_size(st, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray8>") {
    Log("Testing RepresentativeTable - SchaeferTable<BitArray8>");
    auto st = SchaeferTable<BitArray8>(4, 2, 3);
    test_invariants(st, irrep4_0);
    check_size(st, irrep4_0, 3);
  }

  // =====================================================================
  // Cross-check: BoundedPartitions vs SchaeferTable enumerate identically
  // =====================================================================

  SECTION("BoundedPartitions vs SchaeferTable cross-check") {
    Log("Testing RepresentativeTable - BoundedPartitions vs SchaeferTable cross-check");
    auto bp = BoundedPartitions<BitArray2>(4, 2, 3);
    auto st = SchaeferTable<BitArray2>(4, 2, 3);

    RepresentativeTable<BoundedPartitions<BitArray2>> tbl_bp(bp, irrep4_0.group(), irrep4_0.characters());
    RepresentativeTable<SchaeferTable<BitArray2>> tbl_st(st, irrep4_0.group(), irrep4_0.characters());
    REQUIRE(tbl_bp.size() == tbl_st.size());

    // Representatives must appear in the same order
    using bit_t = typename BoundedPartitions<BitArray2>::bit_t;
    std::vector<bit_t> reps_bp, reps_st;
    for (auto r : tbl_bp)
      reps_bp.push_back(r);
    for (auto r : tbl_st)
      reps_st.push_back(r);
    REQUIRE(reps_bp == reps_st);
  }

  // =====================================================================
  // Cross-check: cyclic_group(3) on BoundedMultisets(3, bound=2)
  //
  // 8 states. Orbits:
  //   {(0,0,0)}, {(1,0,0),(0,1,0),(0,0,1)},
  //   {(1,1,0),(0,1,1),(1,0,1)}, {(1,1,1)}
  // k=0: 4 reps
  // k=1: {(0,0,0)} and {(1,1,1)} excluded (chi-sum over stabilizer=G is 0)
  //   → 2 reps
  // =====================================================================

  SECTION("BoundedMultisets<BitArray1> cyclic_group(3)") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArray1> cyclic_group(3)");
    auto ms = BoundedMultisets<BitArray1>(3, 2);
    test_invariants(ms, irrep3_0);
    test_invariants(ms, irrep3_1);
    check_size(ms, irrep3_0, 4);
    check_size(ms, irrep3_1, 2);
  }

  SECTION("BoundedMultisets<BitArrayLong1> cyclic_group(3)") {
    Log("Testing RepresentativeTable - BoundedMultisets<BitArrayLong1> cyclic_group(3)");
    auto ms = BoundedMultisets<BitArrayLong1>(3, 2);
    test_invariants(ms, irrep3_0);
    test_invariants(ms, irrep3_1);
    check_size(ms, irrep3_0, 4);
    check_size(ms, irrep3_1, 2);
  }

  // =====================================================================
  // Triangular lattice 12-site space group (72 elements, non-abelian)
  //
  // Irreps tested:
  //   Gamma.C6.A  — real, trivial
  //   Gamma.C6.B  — real, non-trivial (chi in {+1,-1})
  //   K.C3.Ea     — complex (K-point, cube roots of unity)
  // =====================================================================

  SECTION("triangular_12_site_group") {
    Log("Testing RepresentativeTable - triangular_12_site_group");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.toml";
    auto fl = FileToml(lfile);

    auto irrep_trivial = read_representation(fl, "Gamma.C6.A");
    auto irrep_B = read_representation(fl, "Gamma.C6.B");
    auto irrep_Kea = read_representation(fl, "K.C3.Ea");

    // Combinations<uint32_t>(12, 2): 66 states
    {
      auto combos = Combinations<uint32_t>(12, 2);
      test_invariants(combos, irrep_trivial);
      test_invariants(combos, irrep_B);
      test_invariants(combos, irrep_Kea);

      RepresentativeTable<Combinations<uint32_t>> tbl(combos, irrep_trivial.group(), irrep_trivial.characters());
      REQUIRE(tbl.size() > 0);
      REQUIRE(tbl.size() < combos.size());
    }

    // LinTable<uint32_t>(12, 6): 924 half-filling states
    {
      auto lt = LinTable<uint32_t>(12, 6);
      test_invariants(lt, irrep_trivial);
      test_invariants(lt, irrep_B);
      test_invariants(lt, irrep_Kea);

      RepresentativeTable<LinTable<uint32_t>> tbl(lt, irrep_trivial.group(), irrep_trivial.characters());
      REQUIRE(tbl.size() > 0);
      REQUIRE(tbl.size() < lt.size());
    }

    // BoundedPartitions<BitArray2>(12,2,3): 78 states
    {
      auto bp = BoundedPartitions<BitArray2>(12, 2, 3);
      test_invariants(bp, irrep_trivial);
      test_invariants(bp, irrep_B);
      test_invariants(bp, irrep_Kea);
    }

    // SchaeferTable<BitArray2>(12,2,3): same 78 states
    {
      auto st = SchaeferTable<BitArray2>(12, 2, 3);
      test_invariants(st, irrep_trivial);
      test_invariants(st, irrep_B);

      auto bp = BoundedPartitions<BitArray2>(12, 2, 3);
      RepresentativeTable<BoundedPartitions<BitArray2>> tbl_bp(bp,
                                                               irrep_trivial.group(), irrep_trivial.characters());
      RepresentativeTable<SchaeferTable<BitArray2>> tbl_st(st, irrep_trivial.group(), irrep_trivial.characters());
      REQUIRE(tbl_bp.size() == tbl_st.size());
    }
  }
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
