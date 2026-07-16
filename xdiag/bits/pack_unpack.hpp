// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/bits/bitarray.hpp>

namespace xdiag::bits {

// returns the q-ary represenation of number as a bitarray with n_slots elements.
// n_slots is required when bit_t = BitsetDynamic (dynamic storage must be
// pre-sized); ignored for fixed-size types (default -1 = "not needed").
template <typename bit_t, int nbits>
BitArray<bit_t, nbits> unpack(int64_t number, int64_t q,
                               int64_t n_slots = -1);

// inverse of unpack: returns the integer whose q-ary representation is the
// first n slots of array. Defined inline (hot: called from SchaeferTable::index
// in the matrix-vector product); taken by const reference to avoid copying the
// BitArray on every call.
template <typename bit_t, int nbits>
inline int64_t pack(BitArray<bit_t, nbits> const &array, int64_t q, int64_t n) {
  int64_t result = 0;
  int64_t base = 1;
  for (int64_t i = 0; i < n; ++i) {
    result += array.get(i) * base;
    base *= q;
  }
  return result;
}

} // namespace xdiag::bits
