// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "subsets.hpp"
#include <limits>

#include <extern/fmt/format.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
Subsets<bit_t>::Subsets(int64_t n) try : n_(n), size_(math::ipow(2, n)) {
  if (n < 0) {
    XDIAG_THROW("Error constructing Subsets: n<0");
  }
  if (n > (int64_t)std::numeric_limits<bit_t>::digits) {
    XDIAG_THROW(fmt::format("Error constructing Subsets: n ({}) is too large "
                            "for bit_t with {} bits",
                            n, std::numeric_limits<bit_t>::digits));
  }
}
XDIAG_CATCH

template <class bit_t> int64_t Subsets<bit_t>::n() const { return n_; }
template <class bit_t> int64_t Subsets<bit_t>::size() const { return size_; }
template <class bit_t> int64_t Subsets<bit_t>::bitwidth() const { return n_; }

template <class bit_t> SubsetsIterator<bit_t> Subsets<bit_t>::begin() const {
  return SubsetsIterator<bit_t>((int64_t)0);
}

template <class bit_t> SubsetsIterator<bit_t> Subsets<bit_t>::end() const {
  return SubsetsIterator<bit_t>(size_);
}

template <class bit_t> bit_t Subsets<bit_t>::operator[](int64_t idx) const {
  return (bit_t)idx;
}
template <class bit_t> int64_t Subsets<bit_t>::index(bit_t bits) const {
  return (int64_t)bits;
}

template <class bit_t>
bool Subsets<bit_t>::operator==(Subsets<bit_t> const &rhs) const {
  return n_ == rhs.n_;
}

template <class bit_t>
bool Subsets<bit_t>::operator!=(Subsets<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class Subsets<uint32_t>;
template class Subsets<uint64_t>;

template <class bit_t>
SubsetsIterator<bit_t>::SubsetsIterator(int64_t idx) : current_(idx) {}

template <class bit_t>
bool SubsetsIterator<bit_t>::operator==(
    const SubsetsIterator<bit_t> &rhs) const {
  return current_ == rhs.current_;
}

template <class bit_t>
bool SubsetsIterator<bit_t>::operator!=(
    const SubsetsIterator<bit_t> &rhs) const {
  return !operator==(rhs);
}

template <class bit_t>
SubsetsIterator<bit_t> &SubsetsIterator<bit_t>::operator++() {
  ++current_;
  return *this;
}

template <class bit_t>
SubsetsIterator<bit_t> &SubsetsIterator<bit_t>::operator+=(int64_t n) {
  current_ += n;
  return *this;
}

template <class bit_t>
SubsetsIterator<bit_t> SubsetsIterator<bit_t>::operator+(int64_t n) const {
  SubsetsIterator copy = *this;
  copy += n;
  return copy;
}

template <class bit_t> bit_t SubsetsIterator<bit_t>::operator*() const {
  return (bit_t)current_;
}

template class SubsetsIterator<uint32_t>;
template class SubsetsIterator<uint64_t>;

} // namespace xdiag::combinatorics
