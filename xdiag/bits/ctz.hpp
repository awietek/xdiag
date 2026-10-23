// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/bits/popcount.hpp>

namespace xdiag::bits {

// Portable fallback: (x & -x) - 1 has exactly the bits below the lowest set
// bit set, so its popcount is the number of trailing zeros (the full width for
// x == 0).
constexpr int swar_ctz_32(uint32_t x) noexcept {
  return swar_popcount_32((x & (~x + 1)) - 1);
}
constexpr int swar_ctz_64(uint64_t x) noexcept {
  return swar_popcount_64((x & (~x + 1)) - 1);
}

// Number of trailing zeros, i.e. the position of the lowest set bit. The
// builtins are undefined for x == 0. Compiler coverage as in utils/likely.hpp.
#if defined(__GNUC__) || defined(__clang__) || \
    defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
constexpr int ctz(uint32_t x) noexcept { return __builtin_ctz(x); }
constexpr int ctz(uint64_t x) noexcept { return __builtin_ctzll(x); }
#else
constexpr int ctz(uint32_t x) noexcept { return swar_ctz_32(x); }
constexpr int ctz(uint64_t x) noexcept { return swar_ctz_64(x); }
#endif

} // namespace xdiag::bits
