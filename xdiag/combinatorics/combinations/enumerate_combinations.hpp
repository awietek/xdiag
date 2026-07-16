// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/bits/zero_one.hpp>

namespace xdiag::combinatorics {

template <typename bit_t> bit_t next_combination(bit_t v) noexcept{
  // Bit twiddling Hack from
  // http://graphics.stanford.edu/~seander/bithacks.html
  // #NextBitPermutation

  // DONT USE FAST VERSION

  // // Fast version (needs __builtin_ctz(v)) (some problem with 0)
  // bit_t t = v | (v - 1); // t gets v's least significant 0
  // return ((t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1)));

  // // Fast version (needs __builtin_ctz(v)) (some problem with 0)
  // int_t t = v | (v - 1); // t gets v's least significant 0
  // return v == 0 ? ~v :((t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) +
  // 1)));

  // Slow version (should work everywhere)
  const bit_t zero = bits::zero<bit_t>();
  const bit_t one = bits::one<bit_t>();
  bit_t t = (v | (v - one)) + one;
  return v == zero ? ~v : t | ((((t & -t) / (v & -v)) >> 1) - one);
}

// next_combination(v, n): like next_combination(v) for integral bit_t, but for
// Bitset types uses only test/set/reset (O(1) each) instead of multi-word
// arithmetic, which is O(nchunks) per operation. Prefer this overload from
// iterators that already store n.
template <typename bit_t>
inline bit_t next_combination(bit_t v, int64_t n) noexcept {
  if constexpr (std::is_integral_v<bit_t>) {
    (void)n;
    return next_combination(v);
  } else {
    // Locate lowest set bit
    int64_t low = 0;
    while (low < n && !v.test(low))
      ++low;
    if (low >= n)
      return v;
    // Measure run of consecutive 1s from low
    int64_t m = 0;
    while (low + m < n && v.test(low + m))
      ++m;
    // Guard: called on the last combination (run reaches bit n)
    if (low + m >= n)
      return v;
    // Perform the update using only set/reset
    for (int64_t i = low; i < low + m; ++i)
      v.reset(i);
    v.set(low + m);
    for (int64_t i = 0; i < m - 1; ++i)
      v.set(i);
    return v;
  }
}

// nth_combination(n, k, idx, width): the optional width argument sets the
// storage bit width of the returned bit_t and defaults to n. It matters only
// for dynamic Bitset, whose chunk count tracks the bit count: enumerating k
// bits over n positions but storing them in a wider (width >= n) bitset lets
// the tJ compressed-dn states share the ups bit width, so the kernels'
// nsites-wide masks and the dn states have matching chunk counts.
template <class bit_t>
bit_t nth_combination(int64_t n, int64_t k, int64_t idx, int64_t width = -1);
template <class bit_t> int64_t rank_combination(bit_t bits, int64_t n);

// template <typename bit_t>
// bit_t get_nth_pattern(int64_t n, int64_t nsites, int64_t nupspins);

// template <typename bit_t>
// int64_t get_n_for_pattern(bit_t pattern, int64_t nsites, int64_t nupspins);

} // namespace xdiag::combinatorics
