// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bitarray.hpp"

#include <algorithm>
#include <xdiag/bits/bitset.hpp>
#include <extern/fmt/format.hpp>

namespace xdiag::bits {

// get() and set() are defined inline in bitarray.hpp.

template <typename bit_t, int nbits>
BitArray<bit_t, nbits>::BitArray(bit_t raw) noexcept : bits_(raw) {}

template <typename bit_t, int nbits>
bit_t BitArray<bit_t, nbits>::raw() const noexcept {
  return bits_;
}

template <typename bit_t, int nbits>
bool BitArray<bit_t, nbits>::operator==(
    BitArray<bit_t, nbits> const &rhs) const noexcept {
  return bits_ == rhs.bits_;
}

template <typename bit_t, int nbits>
bool BitArray<bit_t, nbits>::operator!=(
    BitArray<bit_t, nbits> const &rhs) const noexcept {
  return !operator==(rhs);
}

template <typename bit_t, int nbits>
bool BitArray<bit_t, nbits>::operator<(
    BitArray<bit_t, nbits> const &rhs) const noexcept {
  return bits_ < rhs.bits_;
}
template <typename bit_t, int nbits>

bool BitArray<bit_t, nbits>::operator<=(
    BitArray<bit_t, nbits> const &rhs) const noexcept {
  return !operator>(rhs);
}
template <typename bit_t, int nbits>

bool BitArray<bit_t, nbits>::operator>(
    BitArray<bit_t, nbits> const &rhs) const noexcept {
  return rhs.operator<(*this);
}
template <typename bit_t, int nbits>

bool BitArray<bit_t, nbits>::operator>=(
    BitArray<bit_t, nbits> const &rhs) const noexcept {
  return !operator<(rhs);
}

template <typename bit_t, int nbits>
std::string to_string(BitArray<bit_t, nbits> const &bits, int64_t size,
                      bool reverse) {
  std::string str;
  if (nbits <= 3) { // numbers from 0 to 7 only
    for (int64_t i = 0; i < size; ++i) {
      str += std::to_string(bits.get(i));
    }
    if (reverse) {
      std::reverse(str.begin(), str.end());
    }
    return str;
  } else { // numbers with at least two digits
    int digits = std::to_string(bitmask<int>(nbits)).size();
    for (int64_t i = 0; i < size; ++i) {
      int bit = reverse ? bits.get(size - 1 - i) : bits.get(i);
      if (i == size - 1) {
        str += fmt::format("{:{}}", bit, digits);
      } else {
        str += fmt::format("{:{}} ", bit, digits);
      }
    }
    return str;
  }
}

template <typename bit_t, int nbits>
std::ostream &operator<<(std::ostream &out,
                         BitArray<bit_t, nbits> const &bits) {
  out << to_string(bits);
  return out;
}

} // namespace xdiag::bits

using namespace xdiag::bits;

// Template instantiations
#define INSTANTIATE_XDIAG_BITS_BITARRAY(BIT_T, NBITS)                          \
  template class xdiag::bits::BitArray<BIT_T, NBITS>;                          \
  template std::string xdiag::bits::to_string(BitArray<BIT_T, NBITS> const &,  \
                                              int64_t, bool);                  \
  template std::ostream &xdiag::bits::operator<<(                              \
      std::ostream &, BitArray<BIT_T, NBITS> const &);

#define INSTANTIATE_XDIAG_BITS_BITARRAY_FOR_NBITS(BIT_T)                       \
  INSTANTIATE_XDIAG_BITS_BITARRAY(BIT_T, 1)                                    \
  INSTANTIATE_XDIAG_BITS_BITARRAY(BIT_T, 2)                                    \
  INSTANTIATE_XDIAG_BITS_BITARRAY(BIT_T, 3)                                    \
  INSTANTIATE_XDIAG_BITS_BITARRAY(BIT_T, 4)                                    \
  INSTANTIATE_XDIAG_BITS_BITARRAY(BIT_T, 5)                                    \
  INSTANTIATE_XDIAG_BITS_BITARRAY(BIT_T, 6)                                    \
  INSTANTIATE_XDIAG_BITS_BITARRAY(BIT_T, 7)                                    \
  INSTANTIATE_XDIAG_BITS_BITARRAY(BIT_T, 8)

// Native integer types
// BEGIN_INSTANTIATION_GROUP(native)
INSTANTIATE_XDIAG_BITS_BITARRAY_FOR_NBITS(uint32_t)
INSTANTIATE_XDIAG_BITS_BITARRAY_FOR_NBITS(uint64_t)
// END_INSTANTIATION_GROUP

// BitArrayLong1..8 = BitArray<BitsetDynamic, 1..8>
// BEGIN_INSTANTIATION_GROUP(bitarray_long)
INSTANTIATE_XDIAG_BITS_BITARRAY_FOR_NBITS(BitsetDynamic)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bitset_static)
INSTANTIATE_XDIAG_BITS_BITARRAY_FOR_NBITS(BitsetStatic1)
INSTANTIATE_XDIAG_BITS_BITARRAY_FOR_NBITS(BitsetStatic2)
INSTANTIATE_XDIAG_BITS_BITARRAY_FOR_NBITS(BitsetStatic4)
INSTANTIATE_XDIAG_BITS_BITARRAY_FOR_NBITS(BitsetStatic8)
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_BITARRAY
