// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/bits/bitset.hpp>

// (un)comment to use intrinsic BMI2 pext/pdef instruction
#define USE_PEXT_PDEP

#if defined(__BMI2__) && defined(USE_PEXT_PDEP)
#if defined(__x86_64__)
#include <immintrin.h>
#elif defined(__arm__)
#include "arm_neon.h"
#endif

#endif

namespace xdiag::bits {

// wrappers around the bmi2 pdep/pext intrinsics
#if defined(__BMI2__) && defined(USE_PEXT_PDEP)
inline uint16_t pdep(uint16_t source, uint16_t mask) noexcept {
  return _pdep_u32((uint32_t)source, (uint32_t)mask);
}
inline uint32_t pdep(uint32_t source, uint32_t mask) noexcept {
  return _pdep_u32(source, mask);
}
inline uint64_t pdep(uint64_t source, uint64_t mask) noexcept {
  return _pdep_u64(source, mask);
}

inline uint16_t pext(uint16_t source, uint16_t mask) noexcept {
  return _pext_u32((uint32_t)source, (uint32_t)mask);
}
inline uint32_t pext(uint32_t source, uint32_t mask) noexcept {
  return _pext_u32(source, mask);
}
inline uint64_t pext(uint64_t source, uint64_t mask) noexcept {
  return _pext_u64(source, mask);
}
#endif

template <typename bit_t>
inline bit_t pdep_fallback(bit_t x, bit_t mask) noexcept {
  bit_t res = 0;
  for (bit_t bb = 1; mask; bb += bb) {
    if (x & bb) {
      res |= mask & (-mask);
    }
    mask &= (mask - 1);
  }
  return res;
}

template <typename bit_t>
inline bit_t pext_fallback(bit_t x, bit_t mask) noexcept {
  bit_t res = 0;
  for (bit_t bb = 1; mask; bb += bb) {
    if (x & mask & -mask) {
      res |= bb;
    }
    mask &= (mask - 1);
  }
  return res;
}

template <typename bit_t> inline bit_t deposit(bit_t x, bit_t mask) noexcept {
#if defined(__BMI2__) && defined(USE_PEXT_PDEP)
  return pdep(x, mask);
// Parallel Bits Deposit implementation w/o BMI2 intrinsics
#else
  return pdep_fallback(x, mask);
#endif
}

template <typename bit_t> inline bit_t extract(bit_t x, bit_t mask) noexcept {
#if defined(__BMI2__) && defined(USE_PEXT_PDEP)
  return pext(x, mask);
// Parallel Bits Extract implementation w/o BMI2 intrinsics
#else
  return pext_fallback(x, mask);
#endif
}

// Overloads for the dynamic Bitset, which has no integer pext/pdep instruction.
// deposit places the low bits of `x` into the set positions of `mask` (result
// width = mask width); extract packs the masked bits of `x` into the low
// positions (result width = popcount(mask)). Both result widths follow from the
// mask, so the 2-argument signature is preserved; these are more specialized
// than the integral templates above, so overload resolution picks them for
// Bitset arguments. The integer pext/pdep templates are never instantiated for
// Bitset (they would not compile).
template <typename chunk_t>
inline Bitset<chunk_t, 0> deposit(Bitset<chunk_t, 0> const &x,
                                  Bitset<chunk_t, 0> const &mask) {
  int64_t width =
      (int64_t)mask.chunks().size() * (int64_t)Bitset<chunk_t, 0>::nchunkbits;
  Bitset<chunk_t, 0> res(width);
  int64_t rank = 0;
  for (int64_t s = 0; s < width; ++s) {
    if (mask.test(s)) {
      if (x.test(rank)) {
        res.set(s);
      }
      ++rank;
    }
  }
  return res;
}

template <typename chunk_t>
inline Bitset<chunk_t, 0> extract(Bitset<chunk_t, 0> const &x,
                                  Bitset<chunk_t, 0> const &mask) {
  int64_t width =
      (int64_t)mask.chunks().size() * (int64_t)Bitset<chunk_t, 0>::nchunkbits;
  Bitset<chunk_t, 0> res(mask.count());
  int64_t rank = 0;
  for (int64_t s = 0; s < width; ++s) {
    if (mask.test(s)) {
      if (x.test(s)) {
        res.set(rank);
      }
      ++rank;
    }
  }
  return res;
}

} // namespace xdiag::bits
