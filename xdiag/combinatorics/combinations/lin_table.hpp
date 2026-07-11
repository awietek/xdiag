// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>

// LinTable<bit_t> is a drop-in replacement for Combinations<bit_t> with O(1)
// index() lookups via the precomputed split-table method of Lin (1990), PRB 42,
// 6561. The n-bit patterns are split into an upper half (n_left = ceil(n/2)
// bits) and a lower half (n_right = n - n_left bits). Two tables are
// precomputed at construction:
//   left_indices_[left]  – cumulative rank offset for all patterns with the
//                          given upper half (size 2^n_left)
//   right_indices_[right] – rank of each right-half pattern within its
//                           popcount class (size 2^n_right)
// so index(bits) = left_indices_[bits >> n_right] + right_indices_[bits &
// right_mask] is a single O(1) table lookup. Iteration and operator[] delegate
// to Combinations<bit_t>.
// Requires: n >= 0, 0 <= k <= n. Instantiated for uint32_t, uint64_t.
//
// Example:
//   LinTable<uint32_t> lt(16, 8);   // C(16,8) = 12870 patterns
//   for (auto bits : lt)
//     int64_t idx = lt.index(bits); // O(1) lookup
namespace xdiag::combinatorics {

template <class bit_tt> class LinTable {
public:
  using bit_t = bit_tt;
  using iterator_t = CombinationsIterator<bit_t>;

  LinTable() = default;
  LinTable(int64_t n, int64_t k);

  int64_t n() const;
  int64_t k() const;
  constexpr int64_t d() const { return 2; } // Local dim per site
  int64_t size() const;
  int64_t bitwidth() const; // number of bits needed to represent state (=n)

  bit_t operator[](int64_t idx) const;
  inline int64_t index(bit_t bits) const {
    return left_indices_[bits >> n_right_] +
           right_indices_[bits::get_range(bits, 0, n_right_)];
  }
  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(LinTable<bit_t> const &rhs) const;
  bool operator!=(LinTable<bit_t> const &rhs) const;

private:
  Combinations<bit_t> combinations_;

  // Lin Table specifics
  int64_t n_left_;
  int64_t n_right_;
  int64_t left_table_size_;
  int64_t right_table_size_;
  std::vector<int64_t> left_indices_;
  std::vector<int64_t> right_indices_;
};

} // namespace xdiag::combinatorics
