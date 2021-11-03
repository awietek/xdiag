#pragma once

#include <cstdint>

// (un)comment to use intrinsic BMI2 pext/pdef instruction
#define HAS_PEXT_PDEP

#if defined(__BMI2__) && defined(HAS_PEXT_PDEP)
#include <immintrin.h>
#endif

namespace hydra::bitops {

// wrappers around the bmi2 pdep/pext intrinsics
#if defined(__BMI2__) && defined(HAS_PEXT_PDEP)
constexpr uint16_t pdep(uint16_t source, uint16_t mask) noexcept {
  return _pdep_u32((uint32_t)source, (uint32_t)mask);
}
constexpr uint32_t pdep(uint32_t source, uint32_t mask) noexcept {
  return _pdep_u32(source, mask);
}
constexpr uint64_t pdep(uint64_t source, uint64_t mask) noexcept {
  return _pdep_u64(source, mask);
}

constexpr uint16_t pext(uint16_t source, uint16_t mask) noexcept {
  return _pext_u32((uint32_t)source, (uint32_t)mask);
}
constexpr uint32_t pext(uint32_t source, uint32_t mask) noexcept {
  return _pext_u32(source, mask);
}
constexpr uint64_t pext(uint64_t source, uint64_t mask) noexcept {
  return _pext_u64(source, mask);
}
#endif

template <typename bit_t>
constexpr bit_t deposit(bit_t x, bit_t mask) noexcept {
#if defined(__BMI2__) && defined(HAS_PEXT_PDEP)
  return pdep(x, mask);
// Parallel Bits Deposit implementation w/o BMI2 intrinsics
#else
  bit_t res = 0;
  for (bit_t bb = 1; mask; bb += bb) {
    if (x & bb) {
      res |= mask & (-mask);
    }
    mask &= (mask - 1);
  }
  return res;
#endif
}

template <typename bit_t>
constexpr bit_t extract(bit_t x, bit_t mask) noexcept {
#if defined(__BMI2__) && defined(HAS_PEXT_PDEP)
  return pext(x, mask);
// Parallel Bits Extract implementation w/o BMI2 intrinsics
#else
  bit_t res = 0;
  for (bit_t bb = 1; mask; bb += bb) {
    if (x & mask & -mask) {
      res |= bb;
    }
    mask &= (mask - 1);
  }
  return res;
#endif
}

} // namespace hydra::bitops
