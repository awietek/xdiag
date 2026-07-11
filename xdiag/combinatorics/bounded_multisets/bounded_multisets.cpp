// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bounded_multisets.hpp"

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/math/log2.hpp>
#include <xdiag/bits/pack_unpack.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::combinatorics {

template <typename bitarray_t>
BoundedMultisets<bitarray_t>::BoundedMultisets(int64_t n, int64_t d) try
    : n_(n), d_(d), size_(math::ipow((int64_t)d, n)) {
  if (n < 0) {
    XDIAG_THROW("Error constructing BoundedMultisets: n<0");
  }
  if (d < 2) {
    XDIAG_THROW("Error constructing BoundedMultisets: d < 2");
  }
  int64_t nbits_for_d = math::ceillog2(d);
  if (nbits_for_d > nbits) {
    XDIAG_THROW(fmt::format("Error constructing BoundedMultisets: d ({}) "
                            "is too large for BitArray type with nbits ({})",
                            d, nbits));
  }
  if (n > bitarray_t::maximum_size) {
    XDIAG_THROW(
        fmt::format("Error constructing BoundedMultisets: size n ({}) "
                    "is too large for BitArray type with maximum_size ({})",
                    n, bitarray_t::maximum_size));
  }
}
XDIAG_CATCH

template <typename bitarray_t> int64_t BoundedMultisets<bitarray_t>::n() const {
  return n_;
}

template <typename bitarray_t> int64_t BoundedMultisets<bitarray_t>::d() const {
  return d_;
}

template <typename bitarray_t>
int64_t BoundedMultisets<bitarray_t>::size() const {
  return size_;
}

template <typename bitarray_t>
int64_t BoundedMultisets<bitarray_t>::bitwidth() const {
  return n_ * nbits;
}

template <typename bitarray_t>
auto BoundedMultisets<bitarray_t>::operator[](int64_t idx) const -> bitarray_t {
  return bits::unpack<raw_t, nbits>(idx, d_, n_);
}

template <typename bitarray_t>
int64_t BoundedMultisets<bitarray_t>::index(bitarray_t seq) const {
  return bits::pack<raw_t, nbits>(seq, d_, n_);
}

template <typename bitarray_t>
BoundedMultisetsIterator<bitarray_t>
BoundedMultisets<bitarray_t>::begin() const {
  return BoundedMultisetsIterator<bitarray_t>(n_, 0, d_);
}

template <typename bitarray_t>
BoundedMultisetsIterator<bitarray_t> BoundedMultisets<bitarray_t>::end() const {
  return BoundedMultisetsIterator<bitarray_t>(n_, size_, d_);
}

template <typename bitarray_t>
bool BoundedMultisets<bitarray_t>::operator==(
    BoundedMultisets<bitarray_t> const &rhs) const {
  return (n_ == rhs.n_) && (d_ == rhs.d_);
}

template <typename bitarray_t>
bool BoundedMultisets<bitarray_t>::operator!=(
    BoundedMultisets<bitarray_t> const &rhs) const {
  return !operator==(rhs);
}

template <typename bitarray_t>
BoundedMultisetsIterator<bitarray_t>::BoundedMultisetsIterator(int64_t n,
                                                               int64_t idx,
                                                               int64_t d)
    : n_(n), idx_(idx), d_(d),
      current_(bits::unpack<bit_t, nbits>(idx, d, n)) {}

template <typename bitarray_t>
bool BoundedMultisetsIterator<bitarray_t>::operator==(
    BoundedMultisetsIterator<bitarray_t> const &rhs) const {
  return idx_ == rhs.idx_;
}

template <typename bitarray_t>
bool BoundedMultisetsIterator<bitarray_t>::operator!=(
    BoundedMultisetsIterator<bitarray_t> const &rhs) const {
  return idx_ != rhs.idx_;
}

template <typename bitarray_t>
BoundedMultisetsIterator<bitarray_t> &
BoundedMultisetsIterator<bitarray_t>::operator++() {
  ++idx_;
  for (int64_t i = 0; i < n_; ++i) {
    int64_t val = current_.get(i) + 1;
    if (val < d_) {
      current_.set(i, val);
      return *this;
    }
    current_.set(i, 0);
  }
  return *this;
}

template <typename bitarray_t>
BoundedMultisetsIterator<bitarray_t> &
BoundedMultisetsIterator<bitarray_t>::operator+=(int64_t n) {
  idx_ += n;
  current_ = bits::unpack<bit_t, nbits>(idx_, d_, n_);
  return *this;
}

template <typename bitarray_t>
BoundedMultisetsIterator<bitarray_t>
BoundedMultisetsIterator<bitarray_t>::operator+(int64_t n) const {
  BoundedMultisetsIterator copy = *this;
  copy += n;
  return copy;
}

template <typename bitarray_t>
auto BoundedMultisetsIterator<bitarray_t>::operator*() const -> bitarray_t {
  return current_;
}

// clang-format off
#define INSTANTIATE_BMS(BIT_T, NBITS)                                          \
  template class BoundedMultisets<bits::BitArray<BIT_T, NBITS>>;              \
  template class BoundedMultisetsIterator<bits::BitArray<BIT_T, NBITS>>;

#define INSTANTIATE_BMS_ALL(BIT_T)                                             \
  INSTANTIATE_BMS(BIT_T, 1)                                                    \
  INSTANTIATE_BMS(BIT_T, 2)                                                    \
  INSTANTIATE_BMS(BIT_T, 3)                                                    \
  INSTANTIATE_BMS(BIT_T, 4)                                                    \
  INSTANTIATE_BMS(BIT_T, 8)

INSTANTIATE_BMS_ALL(uint32_t)
INSTANTIATE_BMS_ALL(uint64_t)
INSTANTIATE_BMS_ALL(bits::BitsetDynamic)
INSTANTIATE_BMS_ALL(bits::BitsetStatic1)
INSTANTIATE_BMS_ALL(bits::BitsetStatic2)
INSTANTIATE_BMS_ALL(bits::BitsetStatic4)
INSTANTIATE_BMS_ALL(bits::BitsetStatic8)

#undef INSTANTIATE_BMS_ALL
#undef INSTANTIATE_BMS
// clang-format on

} // namespace xdiag::combinatorics
