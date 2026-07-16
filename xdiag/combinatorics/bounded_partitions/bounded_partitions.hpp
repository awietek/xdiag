// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::combinatorics {

template <typename bitarray_t> class BoundedPartitionsIterator;

// BoundedPartitions<bitarray_t> enumerates all ordered sequences of length n
// with elements in {0, ..., d-1} whose sum equals total (weak compositions
// of total into n bounded parts). Each sequence is packed into a
// bitarray_t = BitArray<bit_t, nbits>.
// Sequences are produced in reverse-lexicographic order (slot n-1 is the
// primary key, slot 0 the fastest-varying), consistent with BoundedMultisets.
// The count is given by inclusion-exclusion:
//   sum_{k=0}^{floor(total/d)} (-1)^k C(n,k) C(total-k*d+n-1, n-1)
// Requires: n >= 0, total >= 0, d >= 2,
//           ceillog2(d) <= nbits, n <= bitarray_t::maximum_size.
//
// Example:
//   using A = BitArray<uint64_t, 2>;     // 2-bit slots -> d up to 4
//   BoundedPartitions<A> bp(3, 4, 3);   // triples from {0,1,2} summing to 4
//   for (auto seq : bp)
//     use seq.get(0), seq.get(1), seq.get(2);
template <typename bitarray_t> class BoundedPartitions {
public:
  using bit_t = bitarray_t;
  using raw_t = typename bitarray_t::bit_t;
  static constexpr int nbits = bitarray_t::nbits;
  using iterator_t = BoundedPartitionsIterator<bitarray_t>;

  BoundedPartitions() = default;
  BoundedPartitions(int64_t n, int64_t total, int64_t d);

  int64_t n() const;
  int64_t total() const;
  int64_t d() const; // Local Hilbert space dimension per site
  int64_t size() const;
  int64_t bitwidth() const; // number of bits needed to represent state
  bitarray_t operator[](int64_t idx) const; // Sequence at index idx
  int64_t index(bitarray_t seq) const;      // Index of given sequence
  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(BoundedPartitions<bitarray_t> const &rhs) const;
  bool operator!=(BoundedPartitions<bitarray_t> const &rhs) const;

private:
  int64_t n_ = 0;
  int64_t total_ = 0;
  int64_t d_ = 0;
  int64_t size_ = 0;
};

template <typename bitarray_t> class BoundedPartitionsIterator {
public:
  static constexpr int nbits = bitarray_t::nbits;

  BoundedPartitionsIterator() = default;
  BoundedPartitionsIterator(int64_t n, int64_t total, int64_t d, int64_t idx,
                            bitarray_t current);

  inline bool
  operator==(BoundedPartitionsIterator<bitarray_t> const &rhs) const {
    return idx_ == rhs.idx_;
  }
  inline bool
  operator!=(BoundedPartitionsIterator<bitarray_t> const &rhs) const {
    return idx_ != rhs.idx_;
  }
  BoundedPartitionsIterator &operator++();
  BoundedPartitionsIterator &operator+=(int64_t n);
  BoundedPartitionsIterator operator+(int64_t n) const;
  inline bitarray_t operator*() const { return current_; }

private:
  int64_t n_ = 0;
  int64_t total_ = 0;
  int64_t d_ = 0;
  int64_t idx_ = 0;
  bitarray_t current_;
};

} // namespace xdiag::combinatorics
