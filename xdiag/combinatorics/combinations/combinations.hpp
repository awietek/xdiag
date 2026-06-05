// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::combinatorics {

template <typename bit_t> class CombinationsIterator;

// Represents all combinations of k bits set in n-bit integers.
//
// This class provides iteration over all bit patterns with exactly k bits set
// out of n total bit positions. The patterns are generated in lexicographic
// order (treating bit patterns as unsigned integers).
//
// Template parameter bit_t: uint16_t, uint32_t, uint64_t, or Bitset
//
// Example usage:
//   Combinations<uint64_t> combs(5, 2);  // All ways to choose 2 from 5
//   for (auto pattern : combs) {
//     // pattern will be: 0b00011, 0b00101, 0b00110, 0b01001, ...
//   }
//
// Note: For Bitset types, ensure n does not exceed the bit capacity.
template <typename bit_tt> class Combinations {
public:
  using bit_t = bit_tt;
  using iterator_t = CombinationsIterator<bit_t>;

  Combinations() = default;

  // Construct combinations of k bits from n positions.
  // Throws if k > n, k < 0, or n < 0.
  Combinations(int64_t n, int64_t k);

  int64_t n() const;                            // Total number of bit positions
  int64_t k() const;                            // Number of bits set in each pattern
  constexpr int64_t d() const { return 2; }     // Local dim per site
  int64_t size() const;                         // Total number of combinations C(n,k)
  int64_t bitwidth() const;                     // bits needed to represent state (=n)
  bit_t operator[](int64_t idx) const; // Bit pattern at index idx
  int64_t index(bit_t bits) const;     // Index of given bit pattern
  iterator_t begin() const;            // Iterator to first combination
  iterator_t end() const;              // Iterator past last combination

  bool operator==(Combinations<bit_t> const &rhs) const;
  bool operator!=(Combinations<bit_t> const &rhs) const;

private:
  int64_t n_ = 0;
  int64_t k_ = 0;
  int64_t size_ = 0;
};

// Forward iterator over combination bit patterns.
//
// Iterates through all bit patterns with exactly k bits set in n positions,
// generating patterns in lexicographic order. Uses the "next permutation"
// bit-twiddling algorithm for efficient iteration.
template <typename bit_tt> class CombinationsIterator {
public:
  using bit_t = bit_tt;
  CombinationsIterator() = default;
  CombinationsIterator(int64_t n, int64_t k, int64_t idx);

  bool operator==(CombinationsIterator<bit_t> const &rhs) const;
  bool operator!=(CombinationsIterator<bit_t> const &rhs) const;
  CombinationsIterator &operator++();
  CombinationsIterator &operator+=(int64_t n);
  CombinationsIterator operator+(int64_t n) const;
  bit_t operator*() const;

private:
  bit_t current_;
  int64_t idx_;
  int64_t n_;
};

} // namespace xdiag::combinatorics
