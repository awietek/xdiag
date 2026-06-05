// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::combinatorics {

template <typename bitarray_t> class BoundedMultisetsIterator;

// BoundedMultisets<bitarray_t> enumerates all ordered sequences of length n
// with elements drawn from {0, ..., d-1}. Each sequence is packed into a
// bitarray_t = BitArray<bit_t, nbits> using nbits bits per slot (runtime d,
// compile-time packing width). The constructor checks that
// ceillog2(d) <= nbits and n <= bitarray_t::maximum_size.
// Sequences are produced in little-endian base-d order (slot 0 is the
// least significant digit); total count is d^n.
// Requires: n >= 0, d >= 2.
//
// Example:
//   using A = BitArray<uint64_t, 2>;     // 2-bit slots -> d up to 4
//   BoundedMultisets<A> ms(3, 3);        // 27 triples from {0, 1, 2}
//   for (auto seq : ms)
//     use seq.get(0), seq.get(1), seq.get(2);
template <typename bitarray_t> class BoundedMultisets {
public:
  using bit_t = bitarray_t;
  using raw_t = typename bitarray_t::bit_t;
  static constexpr int nbits = bitarray_t::nbits;
  using iterator_t = BoundedMultisetsIterator<bitarray_t>;

  BoundedMultisets() = default;
  BoundedMultisets(int64_t n, int64_t d);

  int64_t n() const;
  int64_t d() const; // Local Hilbert space dimension per site
  int64_t size() const;
  int64_t bitwidth() const;

  bitarray_t operator[](int64_t idx) const; // Sequence at index idx
  int64_t index(bitarray_t seq) const;      // Index of given sequence
  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(BoundedMultisets<bitarray_t> const &rhs) const;
  bool operator!=(BoundedMultisets<bitarray_t> const &rhs) const;

private:
  int64_t n_ = 0;
  int64_t d_ = 0;
  int64_t size_ = 0;
};

template <typename bitarray_t> class BoundedMultisetsIterator {
public:
  using bit_t = typename bitarray_t::bit_t;
  static constexpr int nbits = bitarray_t::nbits;

  BoundedMultisetsIterator() = default;
  BoundedMultisetsIterator(int64_t n, int64_t idx, int64_t d);

  bool operator==(BoundedMultisetsIterator<bitarray_t> const &rhs) const;
  bool operator!=(BoundedMultisetsIterator<bitarray_t> const &rhs) const;
  BoundedMultisetsIterator &operator++();
  BoundedMultisetsIterator &operator+=(int64_t n);
  BoundedMultisetsIterator operator+(int64_t n) const;
  bitarray_t operator*() const;

private:
  int64_t n_ = 0;
  int64_t idx_ = 0;
  int64_t d_ = 0;
  bitarray_t current_;
};

} // namespace xdiag::combinatorics
