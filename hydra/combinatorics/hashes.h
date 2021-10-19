#pragma once

#include <hydra/common.h>

namespace hydra::combinatorics {

// Fowler–Noll–Vo hash function for uint32
const uint32 fnv1_prime_uint32 = 0x01000193;
const uint32 fnv1_offset_uint32 = 0x811c9dc5;
constexpr uint32 fnv1_mask_uint32 = fnv1_prime_uint32 * fnv1_offset_uint32;
inline uint32 hash_fnv1(uint32 bits) { return fnv1_mask_uint32 ^ bits; }
inline uint16 hash_fnv1(uint16 bits) { return (uint16)hash_fnv1((uint32)bits); }

// Fowler–Noll–Vo hash function for uint64
const uint64 fnv1_prime_uint64 = 0x00000100000001B3;
const uint64 fnv1_offset_uint64 = 0xcbf29ce484222645;
constexpr uint64 fnv1_mask_uint64 = fnv1_prime_uint64 * fnv1_offset_uint64;
inline uint64 hash_fnv1(uint64 bits) { return fnv1_mask_uint64 ^ bits; }

} // namespace hydra::combinatorics
