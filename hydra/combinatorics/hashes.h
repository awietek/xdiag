#pragma once

#include <hydra/common.h>

namespace hydra::combinatorics {

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
constexpr uint64_t hash_fnv1(uint64_t bits) {
  return fnv1_mask_uint64_t ^ bits;
}

} // namespace hydra::combinatorics
