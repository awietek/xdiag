// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "isrepresentative.hpp"

#include <cstdint>

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>

namespace xdiag::symmetries {

// determines whether a state is a representative
template <typename bit_t>
bool isrepresentative(bit_t state, SitePermutation const &action) {
  for (int64_t sym = 0; sym < action.size(); ++sym) {
    bit_t tstate = action.apply(sym, state);
    if (tstate < state) {
      return false;
    }
  }
  return true;
}

#define INSTANTIATE_ISREPRESENTATIVE(BIT_TYPE)                                 \
  template bool isrepresentative(BIT_TYPE, SitePermutation const &);

using namespace bits;

INSTANTIATE_ISREPRESENTATIVE(uint16_t);
INSTANTIATE_ISREPRESENTATIVE(uint32_t);
INSTANTIATE_ISREPRESENTATIVE(uint64_t);
INSTANTIATE_ISREPRESENTATIVE(BitsetDynamic);
INSTANTIATE_ISREPRESENTATIVE(BitsetStatic1);
INSTANTIATE_ISREPRESENTATIVE(BitsetStatic2);
INSTANTIATE_ISREPRESENTATIVE(BitsetStatic4);
INSTANTIATE_ISREPRESENTATIVE(BitsetStatic8);

INSTANTIATE_ISREPRESENTATIVE(BitArray1);
INSTANTIATE_ISREPRESENTATIVE(BitArray2);
INSTANTIATE_ISREPRESENTATIVE(BitArray3);
INSTANTIATE_ISREPRESENTATIVE(BitArray4);
INSTANTIATE_ISREPRESENTATIVE(BitArray8);

INSTANTIATE_ISREPRESENTATIVE(BitArrayLong1);
INSTANTIATE_ISREPRESENTATIVE(BitArrayLong2);
INSTANTIATE_ISREPRESENTATIVE(BitArrayLong3);
INSTANTIATE_ISREPRESENTATIVE(BitArrayLong4);
INSTANTIATE_ISREPRESENTATIVE(BitArrayLong8);

#undef INSTANTIATE_ISREPRESENTATIVE

} // namespace xdiag::symmetries
