// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <type_traits>

#include <xdiag/bits/bitset.hpp>

namespace xdiag::bits {

// Helper functions for generic bit operations (zero overhead for native types)

// Construct an all-zero / single-low-bit configuration.
template <typename bit_t> constexpr bit_t zero(int64_t nbits = 64) noexcept {
  if constexpr (std::is_integral<bit_t>::value) {
    return 0;
  } else {
    return bit_t(nbits);
  }
}

template <typename bit_t> constexpr bit_t one(int64_t nbits = 64) noexcept {
  if constexpr (std::is_integral<bit_t>::value) {
    return 1;
  } else {
    bit_t result(nbits);
    result.set(0);
    return result;
  }
}

// Test whether any / no bit is set.
constexpr bool nonzero(uint8_t s) { return (bool)s; }
constexpr bool nonzero(uint16_t s) { return (bool)s; }
constexpr bool nonzero(uint32_t s) { return (bool)s; }
constexpr bool nonzero(uint64_t s) { return (bool)s; }
constexpr bool nonzero(int8_t s) { return (bool)s; }
constexpr bool nonzero(int16_t s) { return (bool)s; }
constexpr bool nonzero(int32_t s) { return (bool)s; }
constexpr bool nonzero(int64_t s) { return (bool)s; }

template <typename chunk_t, int64_t n_chunks>
constexpr bool nonzero(Bitset<chunk_t, n_chunks> const &bits) {
  return bits.any();
}

// All bits zero -- the complement of nonzero, for every bit_t.
template <typename bit_t> constexpr bool iszero(bit_t const &bits) {
  return !nonzero(bits);
}

} // namespace xdiag::bits
