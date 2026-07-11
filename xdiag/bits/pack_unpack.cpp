// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "pack_unpack.hpp"

#include <xdiag/bits/bitset.hpp>

namespace xdiag::bits {

// returns the nbits-ary represenation of number as a bitarray
template <typename bit_t, int nbits>
BitArray<bit_t, nbits> unpack(int64_t number, int64_t q, int64_t n_slots) {
  BitArray<bit_t, nbits> array = make_bitarray<BitArray<bit_t, nbits>>(n_slots);
  int64_t idx = 0;
  while (number) {
    array.set(idx++, number % q);
    number /= q;
  }
  return array;
}

template <typename bit_t, int nbits>
int64_t pack(BitArray<bit_t, nbits> array, int64_t q, int64_t n) {
  int64_t result = 0;
  int64_t base = 1;
  for (int64_t i = 0; i < n; ++i) {
    result += array.get(i) * base;
    base *= q;
  }
  return result;
}

#define INSTANTIATE_UNPACK(bit_t, n)                                           \
  template BitArray<bit_t, n> unpack<bit_t, n>(int64_t, int64_t, int64_t);    \
  template int64_t pack<bit_t, n>(BitArray<bit_t, n>, int64_t, int64_t);
#define INSTANTIATE_UNPACK_ALL(bit_t)                                          \
  INSTANTIATE_UNPACK(bit_t, 1)                                                 \
  INSTANTIATE_UNPACK(bit_t, 2)                                                 \
  INSTANTIATE_UNPACK(bit_t, 3)                                                 \
  INSTANTIATE_UNPACK(bit_t, 4)                                                 \
  INSTANTIATE_UNPACK(bit_t, 5)                                                 \
  INSTANTIATE_UNPACK(bit_t, 6)                                                 \
  INSTANTIATE_UNPACK(bit_t, 7)                                                 \
  INSTANTIATE_UNPACK(bit_t, 8)

INSTANTIATE_UNPACK_ALL(uint32_t)
INSTANTIATE_UNPACK_ALL(uint64_t)
INSTANTIATE_UNPACK_ALL(BitsetDynamic)
INSTANTIATE_UNPACK_ALL(BitsetStatic1)
INSTANTIATE_UNPACK_ALL(BitsetStatic2)
INSTANTIATE_UNPACK_ALL(BitsetStatic4)
INSTANTIATE_UNPACK_ALL(BitsetStatic8)

#undef INSTANTIATE_UNPACK_ALL
#undef INSTANTIATE_UNPACK

} // namespace xdiag::bits
