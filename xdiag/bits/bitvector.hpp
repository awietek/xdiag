// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>
#include <iterator>

#include <xdiag/bits/bitset.hpp>

namespace xdiag::bits {

template <typename bit_t> class BitVectorReference;
template <typename bit_t> class BitVectorIterator;
template <typename bit_t> class BitVectorConstIterator;

template <typename bit_t> class BitVector {
public:
  using value_type = bit_t;
  using reference = BitVectorReference<bit_t>;
  using const_reference = bit_t;
  using iterator = BitVectorIterator<bit_t>;
  using const_iterator = BitVectorConstIterator<bit_t>;

  BitVector() = default;
  BitVector(int64_t size, int64_t nbits);

  // Size information
  inline int64_t size() const noexcept { return size_; }
  inline int64_t nbits() const noexcept { return nbits_; }
  inline bool empty() const noexcept { return size_ == 0; }

  // Element access (noexcept versions don't check bounds)
  inline bit_t operator[](int64_t index) const noexcept {
    return storage_.get_range(index * nbits_, nbits_);
  }
  inline BitVectorReference<bit_t> operator[](int64_t index) noexcept {
    return BitVectorReference<bit_t>(nbits_, index, storage_);
  }

  // Bounds-checked access (may throw)
  bit_t at(int64_t index) const;
  BitVectorReference<bit_t> at(int64_t index);

  inline bit_t front() const noexcept { return operator[](0); }
  inline bit_t back() const noexcept { return operator[](size_ - 1); }

  // Iterators
  iterator begin() noexcept;
  iterator end() noexcept;
  const_iterator begin() const noexcept;
  const_iterator end() const noexcept;
  const_iterator cbegin() const noexcept;
  const_iterator cend() const noexcept;

  bool operator==(BitVector<bit_t> const &rhs) const noexcept;
  bool operator!=(BitVector<bit_t> const &rhs) const noexcept;

private:
  int64_t nbits_ = 0;
  int64_t size_ = 0;
  Bitset<bit_t> storage_; // dynamically sized Bitset
};

template <typename bit_t> class BitVectorReference {
public:
  inline BitVectorReference(int64_t nbits, int64_t index,
                            Bitset<bit_t> &storage)
      : nbits_(nbits), index_(index), storage_(storage) {}
  inline BitVectorReference &operator=(bit_t bits) {
    storage_.set_range(index_ * nbits_, nbits_, bits);
    return *this;
  }
  inline operator bit_t() const {
    return storage_.get_range(index_ * nbits_, nbits_);
  }

private:
  int64_t nbits_;
  int64_t index_;
  Bitset<bit_t> &storage_;
};

template <typename bit_t> class BitVectorConstIterator {
public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = bit_t;
  using difference_type = int64_t;
  using pointer = bit_t const *;
  using reference = bit_t;

  BitVectorConstIterator(BitVector<bit_t> const *vec, int64_t index) noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type n) const noexcept;

  BitVectorConstIterator &operator++() noexcept;
  BitVectorConstIterator operator++(int) noexcept;
  BitVectorConstIterator &operator--() noexcept;
  BitVectorConstIterator operator--(int) noexcept;

  BitVectorConstIterator &operator+=(difference_type n) noexcept;
  BitVectorConstIterator &operator-=(difference_type n) noexcept;

  BitVectorConstIterator operator+(difference_type n) const noexcept;
  BitVectorConstIterator operator-(difference_type n) const noexcept;
  difference_type operator-(BitVectorConstIterator const &other) const noexcept;

  bool operator==(BitVectorConstIterator const &other) const noexcept;
  bool operator!=(BitVectorConstIterator const &other) const noexcept;
  bool operator<(BitVectorConstIterator const &other) const noexcept;
  bool operator<=(BitVectorConstIterator const &other) const noexcept;
  bool operator>(BitVectorConstIterator const &other) const noexcept;
  bool operator>=(BitVectorConstIterator const &other) const noexcept;

private:
  BitVector<bit_t> const *vec_;
  int64_t index_;
};

template <typename bit_t> class BitVectorIterator {
public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = bit_t;
  using difference_type = int64_t;
  using pointer = BitVectorReference<bit_t> *;
  using reference = BitVectorReference<bit_t>;

  BitVectorIterator(BitVector<bit_t> *vec, int64_t index) noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type n) const noexcept;

  BitVectorIterator &operator++() noexcept;
  BitVectorIterator operator++(int) noexcept;
  BitVectorIterator &operator--() noexcept;
  BitVectorIterator operator--(int) noexcept;

  BitVectorIterator &operator+=(difference_type n) noexcept;
  BitVectorIterator &operator-=(difference_type n) noexcept;

  BitVectorIterator operator+(difference_type n) const noexcept;
  BitVectorIterator operator-(difference_type n) const noexcept;
  difference_type operator-(BitVectorIterator const &other) const noexcept;

  bool operator==(BitVectorIterator const &other) const noexcept;
  bool operator!=(BitVectorIterator const &other) const noexcept;
  bool operator<(BitVectorIterator const &other) const noexcept;
  bool operator<=(BitVectorIterator const &other) const noexcept;
  bool operator>(BitVectorIterator const &other) const noexcept;
  bool operator>=(BitVectorIterator const &other) const noexcept;

private:
  BitVector<bit_t> *vec_;
  int64_t index_;
};

// Non-member operator+ for iterator + n
template <typename bit_t>
BitVectorConstIterator<bit_t>
operator+(typename BitVectorConstIterator<bit_t>::difference_type n,
          BitVectorConstIterator<bit_t> const &it) {
  return it + n;
}

template <typename bit_t>
BitVectorIterator<bit_t>
operator+(typename BitVectorIterator<bit_t>::difference_type n,
          BitVectorIterator<bit_t> const &it) {
  return it + n;
}

} // namespace xdiag::bits
