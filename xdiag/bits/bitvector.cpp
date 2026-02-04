// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bitvector.hpp"

#include <limits>
#include <sstream>

namespace xdiag::bits {

template <typename bit_t>
BitVector<bit_t>::BitVector(int64_t size, int64_t nbits) try
    : nbits_(nbits), size_(size), storage_(nbits * size) {
  if (nbits <= 0) {
    XDIAG_THROW("Number of bits must be positive");
  }
  if (size < 0) {
    XDIAG_THROW("Size must be non-negative");
  }
  if (nbits > storage_.nchunkbits()) {
    XDIAG_THROW("Number of bits requested larger than maximal number of bits "
                "for chunk type");
  }
}
XDIAG_CATCH

// Element access with bounds checking
template <typename bit_t>
bit_t BitVector<bit_t>::at(int64_t index) const try {
  if (index < 0 || index >= size_) {
    std::ostringstream msg;
    msg << "BitVector::at: index " << index << " out of range [0, " << size_ << ")";
    XDIAG_THROW(msg.str());
  }
  return operator[](index);
}
XDIAG_CATCH

template <typename bit_t>
BitVectorReference<bit_t> BitVector<bit_t>::at(int64_t index) try {
  if (index < 0 || index >= size_) {
    std::ostringstream msg;
    msg << "BitVector::at: index " << index << " out of range [0, " << size_ << ")";
    XDIAG_THROW(msg.str());
  }
  return operator[](index);
}
XDIAG_CATCH

// Iterators
template <typename bit_t>
typename BitVector<bit_t>::iterator BitVector<bit_t>::begin() noexcept {
  return iterator(this, 0);
}

template <typename bit_t>
typename BitVector<bit_t>::iterator BitVector<bit_t>::end() noexcept {
  return iterator(this, size_);
}

template <typename bit_t>
typename BitVector<bit_t>::const_iterator
BitVector<bit_t>::begin() const noexcept {
  return const_iterator(this, 0);
}

template <typename bit_t>
typename BitVector<bit_t>::const_iterator
BitVector<bit_t>::end() const noexcept {
  return const_iterator(this, size_);
}

template <typename bit_t>
typename BitVector<bit_t>::const_iterator
BitVector<bit_t>::cbegin() const noexcept {
  return const_iterator(this, 0);
}

template <typename bit_t>
typename BitVector<bit_t>::const_iterator
BitVector<bit_t>::cend() const noexcept {
  return const_iterator(this, size_);
}

template <typename bit_t>
bool BitVector<bit_t>::operator==(BitVector<bit_t> const &rhs) const noexcept {
  return (nbits_ == rhs.nbits_) && (size_ == rhs.size_) &&
         (storage_ == rhs.storage_);
}

template <typename bit_t>
bool BitVector<bit_t>::operator!=(BitVector<bit_t> const &rhs) const noexcept {
  return !operator==(rhs);
}

// BitVectorConstIterator implementations
template <typename bit_t>
BitVectorConstIterator<bit_t>::BitVectorConstIterator(
    BitVector<bit_t> const *vec, int64_t index) noexcept
    : vec_(vec), index_(index) {}

template <typename bit_t>
typename BitVectorConstIterator<bit_t>::reference
BitVectorConstIterator<bit_t>::operator*() const noexcept {
  return (*vec_)[index_];
}

template <typename bit_t>
typename BitVectorConstIterator<bit_t>::reference
BitVectorConstIterator<bit_t>::operator[](difference_type n) const noexcept {
  return (*vec_)[index_ + n];
}

template <typename bit_t>
BitVectorConstIterator<bit_t> &
BitVectorConstIterator<bit_t>::operator++() noexcept {
  ++index_;
  return *this;
}

template <typename bit_t>
BitVectorConstIterator<bit_t>
BitVectorConstIterator<bit_t>::operator++(int) noexcept {
  BitVectorConstIterator tmp = *this;
  ++index_;
  return tmp;
}

template <typename bit_t>
BitVectorConstIterator<bit_t> &
BitVectorConstIterator<bit_t>::operator--() noexcept {
  --index_;
  return *this;
}

template <typename bit_t>
BitVectorConstIterator<bit_t>
BitVectorConstIterator<bit_t>::operator--(int) noexcept {
  BitVectorConstIterator tmp = *this;
  --index_;
  return tmp;
}

template <typename bit_t>
BitVectorConstIterator<bit_t> &
BitVectorConstIterator<bit_t>::operator+=(difference_type n) noexcept {
  index_ += n;
  return *this;
}

template <typename bit_t>
BitVectorConstIterator<bit_t> &
BitVectorConstIterator<bit_t>::operator-=(difference_type n) noexcept {
  index_ -= n;
  return *this;
}

template <typename bit_t>
BitVectorConstIterator<bit_t>
BitVectorConstIterator<bit_t>::operator+(difference_type n) const noexcept {
  return BitVectorConstIterator(vec_, index_ + n);
}

template <typename bit_t>
BitVectorConstIterator<bit_t>
BitVectorConstIterator<bit_t>::operator-(difference_type n) const noexcept {
  return BitVectorConstIterator(vec_, index_ - n);
}

template <typename bit_t>
typename BitVectorConstIterator<bit_t>::difference_type
BitVectorConstIterator<bit_t>::operator-(
    BitVectorConstIterator const &other) const noexcept {
  return index_ - other.index_;
}

template <typename bit_t>
bool BitVectorConstIterator<bit_t>::operator==(
    BitVectorConstIterator const &other) const noexcept {
  return index_ == other.index_;
}

template <typename bit_t>
bool BitVectorConstIterator<bit_t>::operator!=(
    BitVectorConstIterator const &other) const noexcept {
  return index_ != other.index_;
}

template <typename bit_t>
bool BitVectorConstIterator<bit_t>::operator<(
    BitVectorConstIterator const &other) const noexcept {
  return index_ < other.index_;
}

template <typename bit_t>
bool BitVectorConstIterator<bit_t>::operator<=(
    BitVectorConstIterator const &other) const noexcept {
  return index_ <= other.index_;
}

template <typename bit_t>
bool BitVectorConstIterator<bit_t>::operator>(
    BitVectorConstIterator const &other) const noexcept {
  return index_ > other.index_;
}

template <typename bit_t>
bool BitVectorConstIterator<bit_t>::operator>=(
    BitVectorConstIterator const &other) const noexcept {
  return index_ >= other.index_;
}

// BitVectorIterator implementations
template <typename bit_t>
BitVectorIterator<bit_t>::BitVectorIterator(BitVector<bit_t> *vec,
                                            int64_t index) noexcept
    : vec_(vec), index_(index) {}

template <typename bit_t>
typename BitVectorIterator<bit_t>::reference
BitVectorIterator<bit_t>::operator*() const noexcept {
  return (*vec_)[index_];
}

template <typename bit_t>
typename BitVectorIterator<bit_t>::reference
BitVectorIterator<bit_t>::operator[](difference_type n) const noexcept {
  return (*vec_)[index_ + n];
}

template <typename bit_t>
BitVectorIterator<bit_t> &BitVectorIterator<bit_t>::operator++() noexcept {
  ++index_;
  return *this;
}

template <typename bit_t>
BitVectorIterator<bit_t> BitVectorIterator<bit_t>::operator++(int) noexcept {
  BitVectorIterator tmp = *this;
  ++index_;
  return tmp;
}

template <typename bit_t>
BitVectorIterator<bit_t> &BitVectorIterator<bit_t>::operator--() noexcept {
  --index_;
  return *this;
}

template <typename bit_t>
BitVectorIterator<bit_t> BitVectorIterator<bit_t>::operator--(int) noexcept {
  BitVectorIterator tmp = *this;
  --index_;
  return tmp;
}

template <typename bit_t>
BitVectorIterator<bit_t> &
BitVectorIterator<bit_t>::operator+=(difference_type n) noexcept {
  index_ += n;
  return *this;
}

template <typename bit_t>
BitVectorIterator<bit_t> &
BitVectorIterator<bit_t>::operator-=(difference_type n) noexcept {
  index_ -= n;
  return *this;
}

template <typename bit_t>
BitVectorIterator<bit_t>
BitVectorIterator<bit_t>::operator+(difference_type n) const noexcept {
  return BitVectorIterator(vec_, index_ + n);
}

template <typename bit_t>
BitVectorIterator<bit_t>
BitVectorIterator<bit_t>::operator-(difference_type n) const noexcept {
  return BitVectorIterator(vec_, index_ - n);
}

template <typename bit_t>
typename BitVectorIterator<bit_t>::difference_type
BitVectorIterator<bit_t>::operator-(
    BitVectorIterator const &other) const noexcept {
  return index_ - other.index_;
}

template <typename bit_t>
bool BitVectorIterator<bit_t>::operator==(
    BitVectorIterator const &other) const noexcept {
  return index_ == other.index_;
}

template <typename bit_t>
bool BitVectorIterator<bit_t>::operator!=(
    BitVectorIterator const &other) const noexcept {
  return index_ != other.index_;
}

template <typename bit_t>
bool BitVectorIterator<bit_t>::operator<(
    BitVectorIterator const &other) const noexcept {
  return index_ < other.index_;
}

template <typename bit_t>
bool BitVectorIterator<bit_t>::operator<=(
    BitVectorIterator const &other) const noexcept {
  return index_ <= other.index_;
}

template <typename bit_t>
bool BitVectorIterator<bit_t>::operator>(
    BitVectorIterator const &other) const noexcept {
  return index_ > other.index_;
}

template <typename bit_t>
bool BitVectorIterator<bit_t>::operator>=(
    BitVectorIterator const &other) const noexcept {
  return index_ >= other.index_;
}

// Explicit template instantiations
template class BitVector<uint8_t>;
template class BitVector<uint16_t>;
template class BitVector<uint32_t>;
template class BitVector<uint64_t>;

template class BitVectorConstIterator<uint8_t>;
template class BitVectorConstIterator<uint16_t>;
template class BitVectorConstIterator<uint32_t>;
template class BitVectorConstIterator<uint64_t>;

template class BitVectorIterator<uint8_t>;
template class BitVectorIterator<uint16_t>;
template class BitVectorIterator<uint32_t>;
template class BitVectorIterator<uint64_t>;

} // namespace xdiag::bits
