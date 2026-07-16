// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bitvector.hpp"

#include <sstream>
#include <xdiag/utils/error.hpp>

namespace xdiag::bits {

// --- BitVector ---

template <typename value_t>
BitVector<value_t>::BitVector(int64_t size, int64_t nbits) try
    : nbits_(nbits), size_(size), storage_(nbits * size) {
  if (nbits <= 0) {
    XDIAG_THROW("Number of bits must be positive");
  }
  if (size < 0) {
    XDIAG_THROW("Size must be non-negative");
  }
}
XDIAG_CATCH

template <typename value_t>
value_t BitVector<value_t>::at(int64_t index) const try {
  if (index < 0 || index >= size_) {
    std::ostringstream msg;
    msg << "BitVector::at: index " << index << " not in range 0.." << size_;
    XDIAG_THROW(msg.str());
  }
  return operator[](index);
}
XDIAG_CATCH

template <typename value_t>
BitVectorReference<value_t> BitVector<value_t>::at(int64_t index) try {
  if (index < 0 || index >= size_) {
    std::ostringstream msg;
    msg << "BitVector::at: index " << index << " not in range 0.." << size_;
    XDIAG_THROW(msg.str());
  }
  return operator[](index);
}
XDIAG_CATCH

template <typename value_t>
typename BitVector<value_t>::iterator BitVector<value_t>::begin() noexcept {
  return iterator(this, 0);
}

template <typename value_t>
typename BitVector<value_t>::iterator BitVector<value_t>::end() noexcept {
  return iterator(this, size_);
}

template <typename value_t>
typename BitVector<value_t>::const_iterator
BitVector<value_t>::begin() const noexcept {
  return const_iterator(this, 0);
}

template <typename value_t>
typename BitVector<value_t>::const_iterator
BitVector<value_t>::end() const noexcept {
  return const_iterator(this, size_);
}

template <typename value_t>
typename BitVector<value_t>::const_iterator
BitVector<value_t>::cbegin() const noexcept {
  return const_iterator(this, 0);
}

template <typename value_t>
typename BitVector<value_t>::const_iterator
BitVector<value_t>::cend() const noexcept {
  return const_iterator(this, size_);
}

template <typename value_t>
bool BitVector<value_t>::operator==(
    BitVector<value_t> const &rhs) const noexcept {
  return (nbits_ == rhs.nbits_) && (size_ == rhs.size_) &&
         (storage_ == rhs.storage_);
}

template <typename value_t>
bool BitVector<value_t>::operator!=(
    BitVector<value_t> const &rhs) const noexcept {
  return !operator==(rhs);
}

// --- BitVectorReference ---

template <typename value_t>
BitVectorReference<value_t>::BitVectorReference(int64_t nbits, int64_t index,
                                                BitsetDynamic &storage)
    : nbits_(nbits), index_(index), storage_(storage) {}

template <typename value_t>
BitVectorReference<value_t> &
BitVectorReference<value_t>::operator=(value_t const &val) {
  set_element<value_t>(storage_, index_ * nbits_, nbits_, val);
  return *this;
}

template <typename value_t>
BitVectorReference<value_t>::operator value_t() const {
  return get_element<value_t>(storage_, index_ * nbits_, nbits_);
}

// --- BitVectorConstIterator ---

template <typename value_t>
BitVectorConstIterator<value_t>::BitVectorConstIterator(
    BitVector<value_t> const *vec, int64_t index) noexcept
    : vec_(vec), index_(index) {}

template <typename value_t>
typename BitVectorConstIterator<value_t>::reference
BitVectorConstIterator<value_t>::operator*() const noexcept {
  return (*vec_)[index_];
}

template <typename value_t>
typename BitVectorConstIterator<value_t>::reference
BitVectorConstIterator<value_t>::operator[](difference_type n) const noexcept {
  return (*vec_)[index_ + n];
}

template <typename value_t>
BitVectorConstIterator<value_t> &
BitVectorConstIterator<value_t>::operator++() noexcept {
  ++index_;
  return *this;
}

template <typename value_t>
BitVectorConstIterator<value_t>
BitVectorConstIterator<value_t>::operator++(int) noexcept {
  BitVectorConstIterator tmp = *this;
  ++index_;
  return tmp;
}

template <typename value_t>
BitVectorConstIterator<value_t> &
BitVectorConstIterator<value_t>::operator--() noexcept {
  --index_;
  return *this;
}

template <typename value_t>
BitVectorConstIterator<value_t>
BitVectorConstIterator<value_t>::operator--(int) noexcept {
  BitVectorConstIterator tmp = *this;
  --index_;
  return tmp;
}

template <typename value_t>
BitVectorConstIterator<value_t> &
BitVectorConstIterator<value_t>::operator+=(difference_type n) noexcept {
  index_ += n;
  return *this;
}

template <typename value_t>
BitVectorConstIterator<value_t> &
BitVectorConstIterator<value_t>::operator-=(difference_type n) noexcept {
  index_ -= n;
  return *this;
}

template <typename value_t>
BitVectorConstIterator<value_t>
BitVectorConstIterator<value_t>::operator+(difference_type n) const noexcept {
  return BitVectorConstIterator(vec_, index_ + n);
}

template <typename value_t>
BitVectorConstIterator<value_t>
BitVectorConstIterator<value_t>::operator-(difference_type n) const noexcept {
  return BitVectorConstIterator(vec_, index_ - n);
}

template <typename value_t>
typename BitVectorConstIterator<value_t>::difference_type
BitVectorConstIterator<value_t>::operator-(
    BitVectorConstIterator const &other) const noexcept {
  return index_ - other.index_;
}

template <typename value_t>
bool BitVectorConstIterator<value_t>::operator==(
    BitVectorConstIterator const &other) const noexcept {
  return index_ == other.index_;
}

template <typename value_t>
bool BitVectorConstIterator<value_t>::operator!=(
    BitVectorConstIterator const &other) const noexcept {
  return index_ != other.index_;
}

template <typename value_t>
bool BitVectorConstIterator<value_t>::operator<(
    BitVectorConstIterator const &other) const noexcept {
  return index_ < other.index_;
}

template <typename value_t>
bool BitVectorConstIterator<value_t>::operator<=(
    BitVectorConstIterator const &other) const noexcept {
  return index_ <= other.index_;
}

template <typename value_t>
bool BitVectorConstIterator<value_t>::operator>(
    BitVectorConstIterator const &other) const noexcept {
  return index_ > other.index_;
}

template <typename value_t>
bool BitVectorConstIterator<value_t>::operator>=(
    BitVectorConstIterator const &other) const noexcept {
  return index_ >= other.index_;
}

// --- BitVectorIterator ---

template <typename value_t>
BitVectorIterator<value_t>::BitVectorIterator(BitVector<value_t> *vec,
                                              int64_t index) noexcept
    : vec_(vec), index_(index) {}

template <typename value_t>
typename BitVectorIterator<value_t>::reference
BitVectorIterator<value_t>::operator*() const noexcept {
  return (*vec_)[index_];
}

template <typename value_t>
typename BitVectorIterator<value_t>::reference
BitVectorIterator<value_t>::operator[](difference_type n) const noexcept {
  return (*vec_)[index_ + n];
}

template <typename value_t>
BitVectorIterator<value_t> &BitVectorIterator<value_t>::operator++() noexcept {
  ++index_;
  return *this;
}

template <typename value_t>
BitVectorIterator<value_t>
BitVectorIterator<value_t>::operator++(int) noexcept {
  BitVectorIterator tmp = *this;
  ++index_;
  return tmp;
}

template <typename value_t>
BitVectorIterator<value_t> &BitVectorIterator<value_t>::operator--() noexcept {
  --index_;
  return *this;
}

template <typename value_t>
BitVectorIterator<value_t>
BitVectorIterator<value_t>::operator--(int) noexcept {
  BitVectorIterator tmp = *this;
  --index_;
  return tmp;
}

template <typename value_t>
BitVectorIterator<value_t> &
BitVectorIterator<value_t>::operator+=(difference_type n) noexcept {
  index_ += n;
  return *this;
}

template <typename value_t>
BitVectorIterator<value_t> &
BitVectorIterator<value_t>::operator-=(difference_type n) noexcept {
  index_ -= n;
  return *this;
}

template <typename value_t>
BitVectorIterator<value_t>
BitVectorIterator<value_t>::operator+(difference_type n) const noexcept {
  return BitVectorIterator(vec_, index_ + n);
}

template <typename value_t>
BitVectorIterator<value_t>
BitVectorIterator<value_t>::operator-(difference_type n) const noexcept {
  return BitVectorIterator(vec_, index_ - n);
}

template <typename value_t>
typename BitVectorIterator<value_t>::difference_type
BitVectorIterator<value_t>::operator-(
    BitVectorIterator const &other) const noexcept {
  return index_ - other.index_;
}

template <typename value_t>
bool BitVectorIterator<value_t>::operator==(
    BitVectorIterator const &other) const noexcept {
  return index_ == other.index_;
}

template <typename value_t>
bool BitVectorIterator<value_t>::operator!=(
    BitVectorIterator const &other) const noexcept {
  return index_ != other.index_;
}

template <typename value_t>
bool BitVectorIterator<value_t>::operator<(
    BitVectorIterator const &other) const noexcept {
  return index_ < other.index_;
}

template <typename value_t>
bool BitVectorIterator<value_t>::operator<=(
    BitVectorIterator const &other) const noexcept {
  return index_ <= other.index_;
}

template <typename value_t>
bool BitVectorIterator<value_t>::operator>(
    BitVectorIterator const &other) const noexcept {
  return index_ > other.index_;
}

template <typename value_t>
bool BitVectorIterator<value_t>::operator>=(
    BitVectorIterator const &other) const noexcept {
  return index_ >= other.index_;
}

} // namespace xdiag::bits

// --- Explicit instantiations ---

#define INSTANTIATE_BITVECTOR(value_t)                                         \
  template class xdiag::bits::BitVector<value_t>;                              \
  template class xdiag::bits::BitVectorReference<value_t>;                     \
  template class xdiag::bits::BitVectorConstIterator<value_t>;                 \
  template class xdiag::bits::BitVectorIterator<value_t>

using namespace xdiag::bits;

// BEGIN_INSTANTIATION_GROUP(native)
INSTANTIATE_BITVECTOR(uint32_t);
INSTANTIATE_BITVECTOR(uint64_t);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bitset)
INSTANTIATE_BITVECTOR(BitsetDynamic);
INSTANTIATE_BITVECTOR(BitsetStatic1);
INSTANTIATE_BITVECTOR(BitsetStatic2);
INSTANTIATE_BITVECTOR(BitsetStatic4);
INSTANTIATE_BITVECTOR(BitsetStatic8);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bitarray)
INSTANTIATE_BITVECTOR(BitArray1);
INSTANTIATE_BITVECTOR(BitArray2);
INSTANTIATE_BITVECTOR(BitArray3);
INSTANTIATE_BITVECTOR(BitArray4);
INSTANTIATE_BITVECTOR(BitArray8);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bitarraylong)
INSTANTIATE_BITVECTOR(BitArrayLong1);
INSTANTIATE_BITVECTOR(BitArrayLong2);
INSTANTIATE_BITVECTOR(BitArrayLong3);
INSTANTIATE_BITVECTOR(BitArrayLong4);
INSTANTIATE_BITVECTOR(BitArrayLong8);
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_BITVECTOR
