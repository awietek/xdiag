// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag::symmetries {

// Applies site permutations from a PermutationGroup to bit states.
// apply(sym, bits) moves the bit at site i to site p[i], where p is the
// sym-th permutation. Supports native integers (uint16/32/64_t) and BitArray.
class SitePermutation {
public:
  XDIAG_API SitePermutation() = default;
  XDIAG_API explicit SitePermutation(PermutationGroup const &group);

  XDIAG_API int64_t size() const;
  XDIAG_API int64_t nsites() const;
  XDIAG_API PermutationGroup const &group() const;
  XDIAG_API bool operator==(SitePermutation const &rhs) const;
  XDIAG_API bool operator!=(SitePermutation const &rhs) const;

  template <typename bit_t>
  XDIAG_API bit_t apply(int64_t sym, bit_t const &bits) const;

  template <typename bit_t, int nbits>
  XDIAG_API bits::BitArray<bit_t, nbits>
  apply(int64_t sym, bits::BitArray<bit_t, nbits> const &bits) const;

private:
  PermutationGroup group_;
  int64_t nsites_;
};

} // namespace xdiag::symmetries
