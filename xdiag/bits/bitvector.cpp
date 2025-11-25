// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bitvector.hpp"

#include <limits>
namespace xdiag::bits {

template <typename bit_t>
BitVector<bit_t>::BitVector(int64_t nbits, int64_t size) try
    : nbits_(nbits), size_(size) {
  int64_t nbitschunk = std::numeric_limits<bit_t>::digits;
  if (nbits > nbitschunk) {
    XDIAG_THROW("Number of bits requested larger than maximal number of bits "
                "for chunk type");
  }
  int64_t nbitstotal = nbits * size;
  int64_t nchunks = std::max(nbitstotal - 1, (int64_t)0) / nbitschunk + 1;
  storage_.resize(nchunks);
}
XDIAG_CATCH

template <typename bit_t>
bool BitVector<bit_t>::operator==(BitVector<bit_t> const &rhs) const {
  return (nbits_ == rhs.nbits_) && (size_ == rhs.size_) &&
         (storage_ == rhs.storage_);
}
template <typename bit_t>
bool BitVector<bit_t>::operator!=(BitVector<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class BitVector<uint8_t>;
template class BitVector<uint16_t>;
template class BitVector<uint32_t>;
template class BitVector<uint64_t>;

} // namespace xdiag::bits
