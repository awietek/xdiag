// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <cstdint>
#include <type_traits>

namespace xdiag::bits {

// Precompiled bit-extract (a.k.a. parallel bits extract / "sheep-and-goats")
// for a FIXED mask, applied many times. This is the use case of compressing
// many `dns` configurations into the same set of non-up sites: the mask
// `not_ups` is constant across the whole inner loop, so the per-mask work is
// precomputed ONCE here and each application is a branchless O(log wordsize)
// (5 steps for 32-bit, 6 for 64-bit), independent of nsites and of popcount.
//
// Unlike bits::extract (which is the hardware `pext` instruction only on
// x86-BMI2 and otherwise an O(nsites) bit-by-bit software loop), this is always
// branchless O(log wordsize) -- the right primitive on platforms without BMI2
// (e.g. ARM / Apple Silicon). The construction follows Hacker's Delight 2nd ed.,
// figure 7-10 ("compress"), split into a one-time precompute (the `mv_` move
// masks) and a cheap apply.
//
// Only the integral bit types (uint32_t / uint64_t) are used; the dynamic
// Bitset path of the symmetric tJ basis falls back to a binary search instead.
// Only the integral bit types (uint32_t / uint64_t) are explicitly instantiated
// and used; MaskCompressor<bit_t> for a non-integral bit_t (the dynamic Bitset)
// is allowed to exist as a type (its constructor / operator() bodies are simply
// never instantiated -- that path falls back to a binary search), so no
// static_assert here.
template <typename bit_t> class MaskCompressor {
public:
  MaskCompressor() = default;

  // Precomputes the move masks for `mask` (defined out-of-line in the .cpp and
  // explicitly instantiated for uint32_t / uint64_t).
  explicit MaskCompressor(bit_t mask);

  // Compress the masked bits of `x` into the low popcount(mask) positions.
  // `x` may carry bits outside the mask; they are discarded (caller decides
  // whether that is an error, e.g. a double occupancy).
  inline bit_t operator()(bit_t x) const {
    x &= mask_;
    for (int i = 0; i < kIters; ++i) {
      bit_t t = x & mv_[i];
      x = (bit_t)((x ^ t) | (t >> (1 << i)));
    }
    return x;
  }

  inline bit_t mask() const { return mask_; }

private:
  static constexpr int kBits = (int)sizeof(bit_t) * 8;
  static constexpr int kIters = (kBits <= 32) ? 5 : 6;

  bit_t mask_{};
  std::array<bit_t, 6> mv_{};
};

} // namespace xdiag::bits
