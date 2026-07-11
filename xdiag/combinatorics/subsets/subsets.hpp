// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <type_traits>

namespace xdiag::combinatorics {

template <typename bit_t> class SubsetsIterator;

// Represents all subsets of n elements (all n-bit integers from 0 to 2^n-1).
//
// This class provides iteration over all 2^n bit patterns of n bit positions,
// i.e., all subsets of an n-element set, in ascending order.
//
// Template parameter bit_t: uint32_t, uint64_t
//
// Example usage:
//   Subsets<uint64_t> subs(3);  // All 8 subsets of a 3-element set
//   for (auto pattern : subs) {
//     // pattern will be: 0b000, 0b001, 0b010, 0b011, ...
//   }
//
// Note: For Bitset types, ensure n does not exceed the bit capacity.
template <typename bit_tt> class Subsets {
public:
  using bit_t = bit_tt;
  using iterator_t = SubsetsIterator<bit_t>;

  Subsets() = default;
  explicit Subsets(int64_t n);

  int64_t n() const;                            // Total number of bit positions
  constexpr int64_t d() const { return 2; }     // Local dim per site
  int64_t size() const;                         // Total number of subsets (2^n)
  int64_t bitwidth() const;                     // bits needed to represent state (=n)
  bit_t operator[](int64_t idx) const; // Bit pattern at index idx
  iterator_t begin() const;            // Iterator to first subset
  iterator_t end() const;              // Iterator past last subset

  int64_t index(bit_t bits) const; // returns the index of given bits

  bool operator==(Subsets<bit_t> const &rhs) const;
  bool operator!=(Subsets<bit_t> const &rhs) const;

private:
  int64_t n_ = 0;
  int64_t size_ = 0;
};

// Forward iterator over subset bit patterns.
//
// Iterates through all bit patterns from 0 to 2^n-1 in ascending order.
template <typename bit_tt> class SubsetsIterator {
public:
  using bit_t = bit_tt;
  SubsetsIterator() = default;
  SubsetsIterator(int64_t idx);

  bool operator==(const SubsetsIterator<bit_t> &rhs) const;
  bool operator!=(const SubsetsIterator<bit_t> &rhs) const;
  SubsetsIterator &operator++();
  SubsetsIterator &operator+=(int64_t n);
  SubsetsIterator operator+(int64_t n) const;
  bit_t operator*() const;

private:
  int64_t current_;
};

// True iff the enumeration is Subsets, i.e. non-number-conserving. Used by the
// number-conservation-agnostic bases (e.g. BasistJ) to branch the np vs no-np
// seam at compile time.
template <typename enumeration_t> struct is_subsets : std::false_type {};
template <typename bit_t>
struct is_subsets<Subsets<bit_t>> : std::true_type {};
template <typename enumeration_t>
inline constexpr bool is_subsets_v = is_subsets<enumeration_t>::value;

} // namespace xdiag::combinatorics
