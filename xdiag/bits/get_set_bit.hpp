// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/bits/bitset.hpp>

// (un)comment to use intrinsic BMI1 bextr instruction
// #define USE_BEXTR

#if defined(__BMI__) && defined(USE_BEXTR)
#include <immintrin.h>
#endif

namespace xdiag::bits {

#if defined(__BMI__) && defined(USE_BEXTR)
constexpr uint32_t bextr(uint16_t src, uint32_t start, uint32_t len) noexcept {
  return _bextr_u32((uint32_t)src, start, len);
}
constexpr uint32_t bextr(uint32_t src, uint32_t start, uint32_t len) noexcept {
  return _bextr_u32(src, start, len);
}
constexpr uint64_t bextr(uint64_t src, uint32_t start, uint32_t len) noexcept {
  return _bextr_u64(src, start, len);
}
#endif

// get_bit
template <class bit_t>
constexpr inline bool get_bit(bit_t x, uint32_t n) noexcept {
#if defined(__BMI__) && defined(USE_BEXTR)
  return bextr(x, n, 1);
#else
  return (x >> n) & (bit_t)1;
#endif
}

template <class chunk_t, int64_t nchunks>
inline bool get_bit(Bitset<chunk_t, nchunks> const &bits, int64_t n) {
  return bits.test(n);
}

// get_bits
template <class bit_t>
constexpr inline bit_t get_bits(bit_t src, uint32_t start,
                                uint32_t len) noexcept {
#if defined(__BMI__) && defined(USE_BEXTR)
  return bextr(src, start, len);
#else
  return (src >> start) & (((bit_t)1 << len) - 1);
#endif
}

// set_bit
template <class bit_t> constexpr void set_bit(bit_t &x, int64_t b) {
  x |= bit_t(1) << b;
}
template <class chunk_t, int64_t nchunks>
inline void set_bit(Bitset<chunk_t, nchunks> &bits, int64_t n) {
  return bits.set(n);
}

template <class bit_t> constexpr void set_bit(bit_t &x, int64_t b, bool val) {
  bit_t mask = bit_t(1) << b;
  x = (x & ~mask) | (bit_t(-val) & mask);
}
template <class chunk_t, int64_t nchunks>
inline void set_bit(Bitset<chunk_t, nchunks> &bits, int64_t n, bool val) {
  return bits.set(n, val);
}

} // namespace xdiag::bits
