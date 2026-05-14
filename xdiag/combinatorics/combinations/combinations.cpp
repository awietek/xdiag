// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "combinations.hpp"

#include <type_traits>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/enumerate_combinations.hpp>
#include <xdiag/math/binomial.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
Combinations<bit_t>::Combinations(int64_t n, int64_t k) try
    : n_(n), k_(k), size_(math::binomial(n, k)) {
  if (n < 0) {
    XDIAG_THROW("Error constructing Combinations: n<0");
  }
  if (k < 0) {
    XDIAG_THROW("Error constructing Combinations: k<0");
  }
  if (k > n) {
    XDIAG_THROW("Error constructing Combinations: k>n");
  }
}
XDIAG_CATCH

template <class bit_t> int64_t Combinations<bit_t>::n() const { return n_; }
template <class bit_t> int64_t Combinations<bit_t>::k() const { return k_; }
template <class bit_t> int64_t Combinations<bit_t>::size() const {
  return size_;
};
template <class bit_t> int64_t Combinations<bit_t>::bitwidth() const {
  return n_;
};

template <class bit_t>
bit_t Combinations<bit_t>::operator[](int64_t idx) const {
  return nth_combination<bit_t>(n_, k_, idx);
}

template <class bit_t> int64_t Combinations<bit_t>::index(bit_t bits) const {
  return rank_combination(bits, n_);
}

template <class bit_t>
CombinationsIterator<bit_t> Combinations<bit_t>::begin() const {
  return CombinationsIterator<bit_t>(n_, k_, 0);
}

template <class bit_t>
CombinationsIterator<bit_t> Combinations<bit_t>::end() const {
  return CombinationsIterator<bit_t>(n_, k_, size_);
}

template <class bit_t>
bool Combinations<bit_t>::operator==(Combinations<bit_t> const &rhs) const {
  return (n_ == rhs.n_) && (k_ == rhs.k_);
}

template <class bit_t>
bool Combinations<bit_t>::operator!=(Combinations<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class Combinations<uint16_t>;
template class Combinations<uint32_t>;
template class Combinations<uint64_t>;

// Bitset instantiations
#define INSTANTIATE_COMBINATIONS(CHUNK_T, NCHUNKS)                             \
  template class Combinations<bits::Bitset<CHUNK_T, NCHUNKS>>;

#define INSTANTIATE_COMBINATIONS_FOR_NCHUNKS(CHUNK_T)                          \
  INSTANTIATE_COMBINATIONS(CHUNK_T, 0)                                         \
  INSTANTIATE_COMBINATIONS(CHUNK_T, 1)                                         \
  INSTANTIATE_COMBINATIONS(CHUNK_T, 2)                                         \
  INSTANTIATE_COMBINATIONS(CHUNK_T, 4)                                         \
  INSTANTIATE_COMBINATIONS(CHUNK_T, 8)

INSTANTIATE_COMBINATIONS_FOR_NCHUNKS(uint8_t)
INSTANTIATE_COMBINATIONS_FOR_NCHUNKS(uint16_t)
INSTANTIATE_COMBINATIONS_FOR_NCHUNKS(uint32_t)
INSTANTIATE_COMBINATIONS_FOR_NCHUNKS(uint64_t)

#undef INSTANTIATE_COMBINATIONS_FOR_NCHUNKS
#undef INSTANTIATE_COMBINATIONS

template <class bit_t>
CombinationsIterator<bit_t>::CombinationsIterator(int64_t n, int64_t k,
                                                  int64_t idx)
    : current_(idx < math::binomial(n, k) ? nth_combination<bit_t>(n, k, idx)
                                          : bit_t{}),
      idx_(idx), n_(n) {}

template <class bit_t>
bool CombinationsIterator<bit_t>::operator==(
    CombinationsIterator<bit_t> const &rhs) const {
  return idx_ == rhs.idx_;
}

template <class bit_t>
bool CombinationsIterator<bit_t>::operator!=(
    CombinationsIterator<bit_t> const &rhs) const {
  return idx_ != rhs.idx_;
}

template <class bit_t>
CombinationsIterator<bit_t> &CombinationsIterator<bit_t>::operator++() {
  current_ = combinatorics::next_combination(current_, n_);
  ++idx_;
  return *this;
}

template <class bit_t>
CombinationsIterator<bit_t> &
CombinationsIterator<bit_t>::operator+=(int64_t n) {
  idx_ += n;
  current_ = nth_combination<bit_t>(n_, bits::popcount(current_), idx_);
  return *this;
}

template <class bit_t>
CombinationsIterator<bit_t>
CombinationsIterator<bit_t>::operator+(int64_t n) const {
  CombinationsIterator copy = *this;
  copy += n;
  return copy;
}

template <class bit_t> bit_t CombinationsIterator<bit_t>::operator*() const {
  return current_;
}

template class CombinationsIterator<uint16_t>;
template class CombinationsIterator<uint32_t>;
template class CombinationsIterator<uint64_t>;

// Bitset instantiations
#define INSTANTIATE_COMBINATIONS_ITERATOR(CHUNK_T, NCHUNKS)                    \
  template class CombinationsIterator<bits::Bitset<CHUNK_T, NCHUNKS>>;

#define INSTANTIATE_COMBINATIONS_ITERATOR_FOR_NCHUNKS(CHUNK_T)                 \
  INSTANTIATE_COMBINATIONS_ITERATOR(CHUNK_T, 0)                                \
  INSTANTIATE_COMBINATIONS_ITERATOR(CHUNK_T, 1)                                \
  INSTANTIATE_COMBINATIONS_ITERATOR(CHUNK_T, 2)                                \
  INSTANTIATE_COMBINATIONS_ITERATOR(CHUNK_T, 4)                                \
  INSTANTIATE_COMBINATIONS_ITERATOR(CHUNK_T, 8)

INSTANTIATE_COMBINATIONS_ITERATOR_FOR_NCHUNKS(uint8_t)
INSTANTIATE_COMBINATIONS_ITERATOR_FOR_NCHUNKS(uint16_t)
INSTANTIATE_COMBINATIONS_ITERATOR_FOR_NCHUNKS(uint32_t)
INSTANTIATE_COMBINATIONS_ITERATOR_FOR_NCHUNKS(uint64_t)

#undef INSTANTIATE_COMBINATIONS_ITERATOR_FOR_NCHUNKS
#undef INSTANTIATE_COMBINATIONS_ITERATOR

} // namespace xdiag::combinatorics
