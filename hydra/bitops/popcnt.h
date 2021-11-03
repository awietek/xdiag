#pragma once

#include <cstdint>

// (un)comment to use intrinsic popcnt instruction
#define HAS_POPCNT

#if defined(HAS_POPCNT)
#include <immintrin.h>
#endif

namespace hydra::bitops {

// Swar popcnts: see
// https://www.chessprogramming.org/Population_Count#SWAR-Popcount See
constexpr int swar_popcnt_32(uint32_t x) {
  x = x - ((x >> 1) & 033333333333) - ((x >> 2) & 011111111111);
  x = (x + (x >> 3)) & 030707070707;
  return x % 63; /* casting out 63 */
}

constexpr uint64_t swar_popcnt_64(uint64_t x) {
  const uint64_t m_1 = 0x5555555555555555LLU;
  const uint64_t m_2 = 0x3333333333333333LLU;
  const uint64_t m_4 = 0x0f0f0f0f0f0f0f0fLLU;
  x = x - ((x >> 1) & m_1);
  x = (x & m_2) + ((x >> 2) & m_2);
  x = (x + (x >> 4)) & m_4;
  return (x * 0x0101010101010101LLU) >> 56;
}

#if defined(HAS_POPCNT)
constexpr int popcnt(int x) noexcept { return __builtin_popcount(x); }
constexpr int popcnt(uint16_t x) noexcept { return __builtin_popcount(x); }
constexpr int popcnt(uint32_t x) noexcept { return __builtin_popcount(x); }
constexpr int popcnt(uint64_t x) noexcept { return __builtin_popcountll(x); }
#else
constexpr int popcnt(int x) { return swar_popcnt_32((uint32_t)x); }
constexpr int popcnt(uint16_t x) { return swar_popcnt_32((uint32_t)x); }
constexpr int popcnt(uint32_t x) { return swar_popcnt_32(x); }
constexpr int popcnt(uint64_t x) { return swar_popcnt_64(x); }
#endif

} // namespace hydra::bitops
