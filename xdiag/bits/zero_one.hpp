// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace xdiag::bits {

// Helper functions for generic bit operations (zero overhead for native types)
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
} // namespace xdiag::bits
