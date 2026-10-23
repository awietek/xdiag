// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/ctz.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag::symmetries {

// Applies site permutations from a PermutationGroup to bit states.
// apply(sym, bits) moves the bit at site i to site p[i], where p is the
// sym-th permutation. Supports native integers (uint32/64_t), Bitset, and
// BitArray.
class SitePermutation {
public:
  SitePermutation() = default;
  explicit SitePermutation(PermutationGroup const &group);

  int64_t size() const;
  int64_t nsites() const;
  PermutationGroup const &group() const;
  bool operator==(SitePermutation const &rhs) const;
  bool operator!=(SitePermutation const &rhs) const;

  inline uint64_t apply(int64_t sym, uint64_t bits) const {
    int64_t const *permutation = group_.ptr(sym);
    uint64_t bitsr = 0;
    while (bits) { // only the occupied sites move: cost ~ popcount(bits)
      bitsr |= uint64_t(1) << permutation[bits::ctz(bits)];
      bits &= bits - 1; // clear the lowest set bit
    }
    return bitsr;
  }

  inline uint32_t apply(int64_t sym, uint32_t bits) const {
    return (uint32_t)apply(sym, (uint64_t)bits);
  }

  template <typename bit_t> bit_t apply(int64_t sym, bit_t const &bits) const;

private:
  PermutationGroup group_;
  int64_t nsites_ = 0;
};

} // namespace xdiag::symmetries
