// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "enumerate_combinations.hpp"

#include <type_traits>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/math/binomial.hpp>

namespace xdiag::combinatorics {

template <typename bit_t> bit_t next_combination(bit_t v) noexcept {
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

// next_combination(v, n): efficient specialization for Bitset types.
// For integral bit_t delegates to next_combination(v).
// For Bitset bit_t replaces multi-word arithmetic with O(1) bit accesses:
//   find the lowest set bit (low), measure the run of consecutive 1s from it
//   (length m), then: clear [low, low+m), set bit low+m, set bits [0, m-1).
// This avoids negation, addition, and division across nchunks words.
template <typename bit_t> bit_t next_combination(bit_t v, int64_t n) noexcept {
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

// rank_combination: colex rank of the k-subset encoded in bits (bits b_0 < b_1
// < ... < b_{k-1} set) = sum_{j=0}^{k-1} binom(b_j, j+1).
template <class bit_t> int64_t rank_combination(bit_t bits, int64_t n) {
  int64_t result = 0, j = 0;
  for (int64_t b = 0; b < n; ++b) {
    if (bits::get(bits, b)) {
      result += math::binomial(b, j + 1);
      ++j;
    }
  }
  return result;
}

// nth_combination: inverse of rank_combination. Greedy descent: for each slot
// j from k down to 1, find the largest b with binom(b,j) <= remaining.
template <class bit_t>
bit_t nth_combination(int64_t n, int64_t k, int64_t idx, int64_t width) {
  // For dynamic Bitset (nchunks==0), default construction yields empty storage;
  // must pass a capacity to allocate enough chunks. The capacity is `width` (>=
  // n), defaulting to n; a wider capacity lets callers (e.g. the tJ compressed-dn
  // fiber) store the k-of-n pattern in a bitset matching a larger bit width. For
  // integers and static Bitsets, the Bitset(width) constructor is equivalent to
  // the default constructor.
  if (width < 0)
    width = n;
  bit_t result;
  if constexpr (std::is_integral_v<bit_t>)
    result = bit_t{};
  else
    result = bit_t(width);
  int64_t remaining = idx;
  for (int64_t j = k; j >= 1; --j) {
    int64_t b = j - 1;
    while (math::binomial(b + 1, j) <= remaining)
      ++b;
    bits::set(result, b);
    remaining -= math::binomial(b, j);
  }
  return result;
}

// template <typename bit_t>
// bit_t get_nth_pattern(int64_t n, int64_t nsites, int64_t nupspins) {
//   bit_t state = bit_zero<bit_t>(nsites + 1);
//   int64_t counter = n;
//   for (int64_t n_varying_bits = nupspins - 1; n_varying_bits >= 0;
//        --n_varying_bits) {
//     int64_t n_combinations = 0;
//     for (int64_t n_allowed_pos = n_varying_bits; n_allowed_pos <= nsites;
//          ++n_allowed_pos) {
//       n_combinations += binomial(n_allowed_pos, n_varying_bits);

//       if (n_combinations > counter) {
//         counter -= n_combinations - binomial(n_allowed_pos, n_varying_bits);
//         state |= (bit_one<bit_t>(nsites + 1) << n_allowed_pos);
//         break;
//       }
//     }
//   }
//   return state;
// }

// template <typename bit_t>
// int64_t get_n_for_pattern(bit_t pattern, int64_t nsites, int64_t nupspins) {
//   int64_t n = 0;
//   bit_t workpattern = pattern;
//   for (int64_t n_varying_bits = nupspins - 1; n_varying_bits >= 0;
//        --n_varying_bits) {
//     for (int64_t i = 0; i <= nsites; ++i) {
//       // MSB is at 2^i
//       if ((bit_one<bit_t>() << (i + 1)) > workpattern) {
//         n += binomial(i, n_varying_bits + 1);
//         workpattern ^= (bit_one<bit_t>() << i);
//         break;
//       }
//     }
//   }
//   return n;
// }

// clang-format off
#define INSTANTIATE_EC(BIT_T)                                                  \
  template BIT_T next_combination<BIT_T>(BIT_T) noexcept;                     \
  template BIT_T next_combination<BIT_T>(BIT_T, int64_t) noexcept;            \
  template BIT_T nth_combination<BIT_T>(int64_t, int64_t, int64_t, int64_t);  \
  template int64_t rank_combination<BIT_T>(BIT_T, int64_t);

INSTANTIATE_EC(uint32_t)
INSTANTIATE_EC(uint64_t)

#define INSTANTIATE_EC_BITSET(CHUNK_T, NCHUNKS)                                \
  template bits::Bitset<CHUNK_T, NCHUNKS>                                      \
  next_combination<bits::Bitset<CHUNK_T, NCHUNKS>>(                            \
      bits::Bitset<CHUNK_T, NCHUNKS>) noexcept;                                \
  template bits::Bitset<CHUNK_T, NCHUNKS>                                      \
  next_combination<bits::Bitset<CHUNK_T, NCHUNKS>>(                            \
      bits::Bitset<CHUNK_T, NCHUNKS>, int64_t) noexcept;                       \
  template bits::Bitset<CHUNK_T, NCHUNKS>                                      \
  nth_combination<bits::Bitset<CHUNK_T, NCHUNKS>>(int64_t, int64_t, int64_t,  \
                                                  int64_t);                    \
  template int64_t                                                              \
  rank_combination<bits::Bitset<CHUNK_T, NCHUNKS>>(                            \
      bits::Bitset<CHUNK_T, NCHUNKS>, int64_t);

#define INSTANTIATE_EC_BITSET_ALL(CHUNK_T)                                     \
  INSTANTIATE_EC_BITSET(CHUNK_T, 0)                                            \
  INSTANTIATE_EC_BITSET(CHUNK_T, 1)                                            \
  INSTANTIATE_EC_BITSET(CHUNK_T, 2)                                            \
  INSTANTIATE_EC_BITSET(CHUNK_T, 4)                                            \
  INSTANTIATE_EC_BITSET(CHUNK_T, 8)

INSTANTIATE_EC_BITSET_ALL(uint32_t)
INSTANTIATE_EC_BITSET_ALL(uint64_t)

#undef INSTANTIATE_EC_BITSET_ALL
#undef INSTANTIATE_EC_BITSET
#undef INSTANTIATE_EC
// clang-format on

} // namespace xdiag::combinatorics
