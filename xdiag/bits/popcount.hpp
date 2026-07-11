// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/bits/bitset.hpp>

// (un)comment to use intrinsic popcount instruction
#define USE_POPCOUNT

#if defined(USE_POPCOUNT)

#if defined(__x86_64__)
#include <immintrin.h>
#elif defined(__arm__)
#include "arm_neon.h"
#endif

#endif

namespace xdiag::bits {

// Swar popcounts: see
// https://www.chessprogramming.org/Population_Count#SWAR-Popcount See
constexpr int swar_popcount_32(uint32_t x) {
  x = x - ((x >> 1) & 033333333333) - ((x >> 2) & 011111111111);
  x = (x + (x >> 3)) & 030707070707;
  return x % 63; /* casting out 63 */
}

constexpr int swar_popcount_64(uint64_t x) {
  const uint64_t m_1 = 0x5555555555555555LLU;
  const uint64_t m_2 = 0x3333333333333333LLU;
  const uint64_t m_4 = 0x0f0f0f0f0f0f0f0fLLU;
  x = x - ((x >> 1) & m_1);
  x = (x & m_2) + ((x >> 2) & m_2);
  x = (x + (x >> 4)) & m_4;
  return (x * 0x0101010101010101LLU) >> 56;
}
constexpr int swar_popcount(uint32_t x) noexcept { return swar_popcount_32(x); }
constexpr int swar_popcount(uint64_t x) noexcept { return swar_popcount_64(x); }

#if defined(USE_POPCOUNT)
constexpr int popcount(int x) noexcept { return __builtin_popcount(x); }
constexpr int popcount(uint32_t x) noexcept { return __builtin_popcount(x); }
constexpr int popcount(uint64_t x) noexcept { return __builtin_popcountll(x); }
#else
constexpr int popcount(int x) noexcept { return swar_popcount_32((uint32_t)x); }
constexpr int popcount(uint32_t x) noexcept { return swar_popcount_32(x); }
constexpr int popcount(uint64_t x) noexcept { return swar_popcount_64(x); }
#endif

// Popcounts for Bitsets
template <typename chunk_t, int64_t nchunks>
constexpr int popcount(Bitset<chunk_t, nchunks> const &bits) noexcept {
  return (int)bits.count();
}

} // namespace xdiag::bits
