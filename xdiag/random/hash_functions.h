#pragma once

#include <xdiag/common.h>

namespace xdiag::random {

// Fowler–Noll–Vo hash function for uint32_t
constexpr uint32_t fnv1_prime_uint32_t = 0x01000193;
constexpr uint32_t fnv1_offset_uint32_t = 0x811c9dc5;
constexpr uint32_t fnv1_mask_uint32_t =
    fnv1_prime_uint32_t * fnv1_offset_uint32_t;
constexpr uint32_t hash_fnv1(uint32_t bits) noexcept {
  return fnv1_mask_uint32_t ^ bits;
}
constexpr uint16_t hash_fnv1(uint16_t bits) noexcept {
  return (uint16_t)hash_fnv1((uint32_t)bits);
}

// Fowler–Noll–Vo hash function for uint64_t
constexpr uint64_t fnv1_prime_uint64_t = 0x00000100000001B3;
constexpr uint64_t fnv1_offset_uint64_t = 0xcbf29ce484222645;
constexpr uint64_t fnv1_mask_uint64_t =
    fnv1_prime_uint64_t * fnv1_offset_uint64_t;
constexpr uint64_t hash_fnv1(uint64_t bits) noexcept {
  return fnv1_mask_uint64_t ^ bits;
}

inline uint64_t hash_div3(uint64_t bits) noexcept {
  uint64_t A = 0;
  uint64_t B = 0;
  uint64_t C = 0;
  int cnt = 0;
  while (bits) {
    A |= (bits & 1) << cnt;
    B |= (bits & 2) << (cnt - 1);
    C |= (bits & 4) << (cnt - 2);
    bits >>= 3;
    ++cnt;
  }
  uint64_t num = A * 1357911 + B * 1197531 + C * 2739651;
  return (123456789 * num + 987654321) % 3000000019ULL;
}

inline uint32_t hash_div3(uint32_t bits) noexcept {
  uint32_t A = 0;
  uint32_t B = 0;
  uint32_t C = 0;
  int cnt = 0;
  while (bits) {
    A |= (bits & 1) << cnt;
    B |= (bits & 2) << (cnt - 1);
    C |= (bits & 4) << (cnt - 2);
    bits >>= 3;
    ++cnt;
  }
  uint32_t num = A * 1357911 + B * 1197531 + C * 2739651;
  return (123456789 * num + 987654321) % 3000000019ULL;
}

inline uint16_t hash_div3(uint16_t bits) noexcept {
  return (uint16_t)hash_div3((uint32_t)bits);
}

// Taken from boost::hash_combine
constexpr uint32_t hash_combine(uint32_t h1, uint32_t h2) {
  h1 ^= h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);
  return h1;
}

constexpr uint64_t hash_combine(uint64_t h1, uint64_t h2) {
  h1 ^= h2 + 0x517cc1b727220a95 + (h1 << 6) + (h1 >> 2);
  return h1;
}

} // namespace xdiag::random
