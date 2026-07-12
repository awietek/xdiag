// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::random {

// Fowler–Noll–Vo hash function for uint64_t
constexpr uint64_t fnv1_prime_uint64_t = 0x00000100000001B3;
constexpr uint64_t fnv1_offset_uint64_t = 0xcbf29ce484222645;
constexpr uint64_t fnv1_mask_uint64_t =
    fnv1_prime_uint64_t * fnv1_offset_uint64_t;
constexpr uint64_t hash_fnv1(uint64_t bits) noexcept {
  return fnv1_mask_uint64_t ^ bits;
}

constexpr uint64_t hash_combine(uint64_t h1, uint64_t h2) {
  h1 ^= h2 + 0x517cc1b727220a95 + (h1 << 6) + (h1 >> 2);
  return h1;
}

// splitmix64 finalizer: a bijective avalanche step that spreads the entropy of
// a hash over all 64 bits. Applied to the public hashes so that low-entropy
// inputs (small charges, nsites, ...) do not leave the high bits nearly
// constant (which made block / irrep IDs share a long common prefix).
constexpr uint64_t hash_finalize(uint64_t h) noexcept {
  h ^= h >> 30;
  h *= 0xbf58476d1ce4e5b9ULL;
  h ^= h >> 27;
  h *= 0x94d049bb133111ebULL;
  h ^= h >> 31;
  return h;
}

} // namespace xdiag::random
