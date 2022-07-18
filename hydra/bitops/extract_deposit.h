#pragma once

#include <cstdint>

// (un)comment to use intrinsic BMI2 pext/pdef instruction
#define USE_PEXT_PDEP

#if defined(__BMI2__) && defined(USE_PEXT_PDEP)

#ifdef __x86_64__
#include <immintrin.h>
#elifdef __arm__
#include "arm_neon.h"
#endif

#endif

namespace hydra::bitops {

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

template <typename bit_t>
inline bit_t deposit(bit_t x, bit_t mask) noexcept {
#if defined(__BMI2__) && defined(USE_PEXT_PDEP)
  return pdep(x, mask);
// Parallel Bits Deposit implementation w/o BMI2 intrinsics
#else
  return pdep_fallback(x, mask);
#endif
}

template <typename bit_t>
inline bit_t extract(bit_t x, bit_t mask) noexcept {
#if defined(__BMI2__) && defined(USE_PEXT_PDEP)
  return pext(x, mask);
// Parallel Bits Extract implementation w/o BMI2 intrinsics
#else
  return pext_fallback(x, mask);
#endif
}

} // namespace hydra::bitops
