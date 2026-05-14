// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "site_permutation.hpp"

#include <type_traits>

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/get_set.hpp>

namespace xdiag::symmetries {

SitePermutation::SitePermutation(PermutationGroup const &group)
    : group_(group), nsites_(group.nsites()) {}

int64_t SitePermutation::size() const { return group_.size(); }
int64_t SitePermutation::nsites() const { return group_.nsites(); }
PermutationGroup const &SitePermutation::group() const { return group_; }

bool SitePermutation::operator==(SitePermutation const &rhs) const {
  return group_ == rhs.group_;
}
bool SitePermutation::operator!=(SitePermutation const &rhs) const {
  return !operator==(rhs);
}

template <typename bit_t>
bit_t SitePermutation::apply(int64_t sym, bit_t const &bits) const {
  bit_t bitsr;
  if constexpr (std::is_same_v<bit_t, bits::BitsetDynamic>) {
    bitsr = bits::BitsetDynamic(nsites_);
  } else {
    bitsr = bit_t{};
  }
  const int64_t *permutation = group_.ptr(sym);
  for (int64_t i = 0; i < nsites_; ++i) {
    bits::set(bitsr, permutation[i], bits::get(bits, i));
  }
  return bitsr;
}

template <typename bit_t, int nbits>
bits::BitArray<bit_t, nbits>
SitePermutation::apply(int64_t sym,
                       bits::BitArray<bit_t, nbits> const &bits) const {
  bits::BitArray<bit_t, nbits> bitsr;
  if constexpr (std::is_same_v<bit_t, bits::BitsetDynamic>) {
    bitsr = bits::BitArray<bit_t, nbits>(bit_t(nsites_ * nbits));
  }
  const int64_t *permutation = group_.ptr(sym);
  for (int64_t i = 0; i < nsites_; ++i) {
    bitsr.set(permutation[i], bits.get(i));
  }
  return bitsr;
}
using namespace bits;

#define INSTANTIATE_APPLY(BIT_TYPE)                                            \
  template BIT_TYPE SitePermutation::apply(int64_t, BIT_TYPE const &) const;

INSTANTIATE_APPLY(uint16_t);
INSTANTIATE_APPLY(uint32_t);
INSTANTIATE_APPLY(uint64_t);
INSTANTIATE_APPLY(BitsetDynamic);
INSTANTIATE_APPLY(BitsetStatic1);
INSTANTIATE_APPLY(BitsetStatic2);
INSTANTIATE_APPLY(BitsetStatic4);
INSTANTIATE_APPLY(BitsetStatic8);

#define INSTANTIATE_APPLY_BITARRAY(BIT_TYPE)                                   \
  template BitArray<BIT_TYPE, 1> SitePermutation::apply(                       \
      int64_t, BitArray<BIT_TYPE, 1> const &) const;                           \
  template BitArray<BIT_TYPE, 2> SitePermutation::apply(                       \
      int64_t, BitArray<BIT_TYPE, 2> const &) const;                           \
  template BitArray<BIT_TYPE, 3> SitePermutation::apply(                       \
      int64_t, BitArray<BIT_TYPE, 3> const &) const;                           \
  template BitArray<BIT_TYPE, 4> SitePermutation::apply(                       \
      int64_t, BitArray<BIT_TYPE, 4> const &) const;                           \
  template BitArray<BIT_TYPE, 5> SitePermutation::apply(                       \
      int64_t, BitArray<BIT_TYPE, 5> const &) const;                           \
  template BitArray<BIT_TYPE, 6> SitePermutation::apply(                       \
      int64_t, BitArray<BIT_TYPE, 6> const &) const;                           \
  template BitArray<BIT_TYPE, 7> SitePermutation::apply(                       \
      int64_t, BitArray<BIT_TYPE, 7> const &) const;                           \
  template BitArray<BIT_TYPE, 8> SitePermutation::apply(                       \
      int64_t, BitArray<BIT_TYPE, 8> const &) const;

INSTANTIATE_APPLY_BITARRAY(uint16_t);
INSTANTIATE_APPLY_BITARRAY(uint32_t);
INSTANTIATE_APPLY_BITARRAY(uint64_t);
INSTANTIATE_APPLY_BITARRAY(BitsetDynamic);
INSTANTIATE_APPLY_BITARRAY(BitsetStatic1);
INSTANTIATE_APPLY_BITARRAY(BitsetStatic2);
INSTANTIATE_APPLY_BITARRAY(BitsetStatic4);
INSTANTIATE_APPLY_BITARRAY(BitsetStatic8);

// #undef INSTANTIATE_APPLY;

} // namespace xdiag::symmetries
