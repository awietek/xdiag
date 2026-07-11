// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "to_string.hpp"
#include <xdiag/bits/bitarray.hpp>

namespace xdiag::bits {

template <typename bit_t>
std::string to_string(bit_t bits, int64_t size, bool reverse) {
  auto ba = BitArray<bit_t, 1>();
  int64_t i = 0;
  while (bits) {
    ba.set(i++, bits & 1);
    bits >>= 1;
  }
  return to_string(ba, size, reverse);
}

template std::string to_string(uint32_t, int64_t, bool);
template std::string to_string(uint64_t, int64_t, bool);

} // namespace xdiag::bits
