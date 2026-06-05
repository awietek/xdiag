// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

#include <xdiag/bits/pack_unpack.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>

// SchaeferTable<bitarray_t> is a drop-in replacement for
// BoundedPartitions<bitarray_t> with fast index() lookups via the precomputed
// split-table method described in Schaefer & Luitz, SciPost Phys. Codebases 48
// (2025), implemented in the DanceQ library.
//
// The n slots are split into a fast subsystem (slots 0..n_fast-1,
// n_fast = ceil(n/2)) and a slow subsystem (slots n_fast..n-1,
// n_slow = floor(n/2)). Two lookup tables are precomputed at construction:
//
//   fast_rank_table_[fast_pack]
//     – local rank of (a[0],...,a[n_fast-1]) within BoundedPartitions(n_fast,
//       k_fast, d), where k_fast = sum of fast-slot values (determined
//       uniquely by fast_pack).
//
//   slow_offset_table_[slow_pack]
//     – cumulative count of all global states whose slow configuration
//       (a[n_fast],...,a[n-1]) precedes the given slow state in rlex order.
//       Built by iterating slow_pack = 0,1,...,d^n_slow-1 (ascending
//       = rlex order) and accumulating count_BP(n_fast, total-slow_sum, d).
//
// Index formula (O(n_slow) to compute slow_pack + O(1) table lookup):
//   index = slow_offset_table_[slow_pack] + fast_rank_table_[fast_pack]
//
// Iteration and operator[] delegate to BoundedPartitions<bitarray_t>.
// Requires: same constraints as BoundedPartitions.
// Instantiated for the same BitArray types as BoundedPartitions.
//
// Example:
//   using A = BitArray<uint64_t, 2>;
//   SchaeferTable<A> st(20, 10, 3);    // n=20, total=10, d=3
//   for (auto seq : st)
//     int64_t idx = st.index(seq);     // fast table lookup

namespace xdiag::combinatorics {

template <typename bitarray_t> class SchaeferTable {
public:
  using bit_t = bitarray_t;
  using raw_t = typename bitarray_t::bit_t;
  static constexpr int nbits = bitarray_t::nbits;
  using iterator_t = BoundedPartitionsIterator<bitarray_t>;

  SchaeferTable() = default;
  SchaeferTable(int64_t n, int64_t total, int64_t d);

  int64_t n() const;
  int64_t total() const;
  int64_t d() const;
  int64_t size() const;
  int64_t bitwidth() const;

  bitarray_t operator[](int64_t idx) const;

  inline int64_t index(bitarray_t seq) const {
    // Fast pack: encode slots 0..n_fast_-1 as a base-d_ integer
    int64_t fast_pack = bits::pack(seq, d_, n_fast_);

    // Slow pack: encode slots n_fast_..n-1 as a base-d_ integer
    // (ascending = rlex order, since a[n_fast_] is the lowest digit)
    int64_t slow_pack = 0, base = 1;
    for (int64_t i = 0; i < n_slow_; ++i) {
      slow_pack += seq.get(n_fast_ + i) * base;
      base *= d_;
    }

    return slow_offset_table_[slow_pack] + fast_rank_table_[fast_pack];
  }

  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(SchaeferTable<bitarray_t> const &rhs) const;
  bool operator!=(SchaeferTable<bitarray_t> const &rhs) const;

private:
  BoundedPartitions<bitarray_t> bounded_partitions_;

  int64_t n_fast_ = 0; // ceil(n/2): number of fast (low-index) slots
  int64_t n_slow_ = 0; // floor(n/2): number of slow (high-index) slots
  int64_t d_ = 0;  // cached for use in inline index()

  // [fast_pack] -> local rank within BoundedPartitions(n_fast_, k_fast, d)
  // k_fast is the digit sum of fast_pack in base d_
  std::vector<int64_t> fast_rank_table_;

  // [slow_pack] -> global index offset for all states with this exact slow
  // configuration; built by iterating slow_pack in ascending order (= rlex)
  std::vector<int64_t> slow_offset_table_;
};

} // namespace xdiag::combinatorics
