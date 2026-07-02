// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <type_traits>

#include <xdiag/bits/zero_one.hpp>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>

namespace xdiag::bits {

template <typename value_t> class BitVectorReference;
template <typename value_t> class BitVectorIterator;
template <typename value_t> class BitVectorConstIterator;

// Dynamic vector of bit-packed values.
//
// BitVector stores a variable number of values (each taking nbits bits)
// packed into a BitsetDynamic (dynamic, always 64-bit chunks).
// Provides a std::vector-like interface with random access and iterators.
//
// Template parameter:
//   value_t: Element type returned by operator[]. Supported types:
//     - uint16_t, uint32_t, uint64_t  (nbits <= type width)
//     - BitsetDynamic, BitsetStatic*   (multi-word values)
//     - BitArray<bit_t, N>             (packed N-bit-per-site arrays;
//                                       pass nbits = digits of bit_t)
//
// Example:
//   BitVector<uint64_t> vec(100, 3);  // 100 elements, 3 bits each (values 0-7)
//   vec[0] = 5;                       // Set first element to 5
//   uint64_t val = vec[0];            // Get first element (returns 5)
//   for (auto x : vec) { ... }        // Iterate over elements

// --- Element access helpers (inline for operator[] performance) ---

template <typename value_t>
inline value_t get_element(BitsetDynamic const &storage, int64_t offset,
                           int64_t nbits) noexcept {
  if constexpr (std::is_integral_v<value_t>) {
    return static_cast<value_t>(storage.get_range(offset, nbits));
  } else if constexpr (detail::is_bitarray_v<value_t>) {
    using raw_t = typename value_t::bit_t;
    return value_t(get_element<raw_t>(storage, offset, nbits));
  } else {
    value_t result = bits::zero<value_t>(nbits);
    for (int64_t b = 0; b < nbits; b += 64) {
      int64_t w = std::min(int64_t(64), nbits - b);
      result.set_range(b, w, storage.get_range(offset + b, w));
    }
    return result;
  }
}

// Atomic OR variant: valid when storage is zero-initialised and each bit
// position is written by exactly one thread (disjoint-orbit guarantee).
template <typename value_t>
inline void set_element_atomic_or(BitsetDynamic &storage, int64_t offset,
                                   int64_t nbits, value_t const &val) noexcept {
  if constexpr (std::is_integral_v<value_t>) {
    storage.set_range_atomic_or(offset, nbits, static_cast<uint64_t>(val));
  } else if constexpr (detail::is_bitarray_v<value_t>) {
    using raw_t = typename value_t::bit_t;
    set_element_atomic_or<raw_t>(storage, offset, nbits, val.raw());
  } else {
    for (int64_t b = 0; b < nbits; b += 64) {
      int64_t w = std::min(int64_t(64), nbits - b);
      storage.set_range_atomic_or(offset + b, w, val.get_range(b, w));
    }
  }
}

template <typename value_t>
inline void set_element(BitsetDynamic &storage, int64_t offset, int64_t nbits,
                        value_t const &val) noexcept {
  if constexpr (std::is_integral_v<value_t>) {
    storage.set_range(offset, nbits, static_cast<uint64_t>(val));
  } else if constexpr (detail::is_bitarray_v<value_t>) {
    using raw_t = typename value_t::bit_t;
    set_element<raw_t>(storage, offset, nbits, val.raw());
  } else {
    for (int64_t b = 0; b < nbits; b += 64) {
      int64_t w = std::min(int64_t(64), nbits - b);
      storage.set_range(offset + b, w, val.get_range(b, w));
    }
  }
}

// --- BitVector ---

template <typename value_t = uint64_t> class BitVector {
public:
  using value_type = value_t;
  using reference = BitVectorReference<value_t>;
  using const_reference = value_t;
  using iterator = BitVectorIterator<value_t>;
  using const_iterator = BitVectorConstIterator<value_t>;

  BitVector() = default;
  BitVector(int64_t size, int64_t nbits);

  inline int64_t size() const noexcept { return size_; }
  inline int64_t nbits() const noexcept { return nbits_; }
  inline bool empty() const noexcept { return size_ == 0; }

  inline value_t operator[](int64_t index) const noexcept {
    return get_element<value_t>(storage_, index * nbits_, nbits_);
  }
  inline BitVectorReference<value_t> operator[](int64_t index) noexcept {
    return BitVectorReference<value_t>(nbits_, index, storage_);
  }

  value_t at(int64_t index) const;
  BitVectorReference<value_t> at(int64_t index);

  // Parallel-safe write: OR the packed bits of val into element index.
  // Requires zero-initialised storage and that each bit position is written
  // by at most one thread (disjoint-orbit guarantee in RepresentativeTable).
  inline void atomic_or_element(int64_t index, value_t const &val) noexcept {
    set_element_atomic_or(storage_, index * nbits_, nbits_, val);
  }

  inline value_t front() const noexcept { return operator[](0); }
  inline value_t back() const noexcept { return operator[](size_ - 1); }

  iterator begin() noexcept;
  iterator end() noexcept;
  const_iterator begin() const noexcept;
  const_iterator end() const noexcept;
  const_iterator cbegin() const noexcept;
  const_iterator cend() const noexcept;

  bool operator==(BitVector<value_t> const &rhs) const noexcept;
  bool operator!=(BitVector<value_t> const &rhs) const noexcept;

private:
  int64_t nbits_ = 0;
  int64_t size_ = 0;
  BitsetDynamic storage_;
};

// --- BitVectorReference ---

template <typename value_t> class BitVectorReference {
public:
  BitVectorReference(int64_t nbits, int64_t index, BitsetDynamic &storage);
  BitVectorReference &operator=(value_t const &val);
  operator value_t() const;

  friend bool operator==(BitVectorReference const &lhs, value_t const &rhs) {
    return static_cast<value_t>(lhs) == rhs;
  }
  friend bool operator!=(BitVectorReference const &lhs, value_t const &rhs) {
    return !(lhs == rhs);
  }
  friend bool operator==(value_t const &lhs, BitVectorReference const &rhs) {
    return lhs == static_cast<value_t>(rhs);
  }
  friend bool operator!=(value_t const &lhs, BitVectorReference const &rhs) {
    return !(lhs == rhs);
  }

private:
  int64_t nbits_;
  int64_t index_;
  BitsetDynamic &storage_;
};

// --- BitVectorConstIterator ---

template <typename value_t> class BitVectorConstIterator {
public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = value_t;
  using difference_type = int64_t;
  using pointer = value_t const *;
  using reference = value_t;

  BitVectorConstIterator(BitVector<value_t> const *vec, int64_t index) noexcept;

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
  BitVector<value_t> const *vec_;
  int64_t index_;
};

// --- BitVectorIterator ---

template <typename value_t> class BitVectorIterator {
public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = value_t;
  using difference_type = int64_t;
  using pointer = BitVectorReference<value_t> *;
  using reference = BitVectorReference<value_t>;

  BitVectorIterator(BitVector<value_t> *vec, int64_t index) noexcept;

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
  BitVector<value_t> *vec_;
  int64_t index_;
};

// Non-member operator+ for n + iterator
template <typename value_t>
BitVectorConstIterator<value_t>
operator+(typename BitVectorConstIterator<value_t>::difference_type n,
          BitVectorConstIterator<value_t> const &it) {
  return it + n;
}

template <typename value_t>
BitVectorIterator<value_t>
operator+(typename BitVectorIterator<value_t>::difference_type n,
          BitVectorIterator<value_t> const &it) {
  return it + n;
}

} // namespace xdiag::bits
