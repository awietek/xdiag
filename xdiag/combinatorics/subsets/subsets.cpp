// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "subsets.hpp"
#include <xdiag/utils/error.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
Subsets<bit_t>::Subsets(int64_t n) try : n_(n), size_((int64_t)1 << n) {
  if (n < 0) {
    XDIAG_THROW("Error constructing Subsets: n<0");
  } else {
    begin_ = SubsetsIterator<bit_t>((int64_t)0);
    end_ = SubsetsIterator<bit_t>((int64_t)1 << n);
  }
}
XDIAG_CATCH

template <class bit_t> int64_t Subsets<bit_t>::n() const { return n_; }
template <class bit_t> int64_t Subsets<bit_t>::size() const { return size_; };
template <class bit_t> SubsetsIterator<bit_t> Subsets<bit_t>::begin() const {
  return begin_;
}
template <class bit_t> SubsetsIterator<bit_t> Subsets<bit_t>::end() const {
  return end_;
}

template class Subsets<uint16_t>;
template class Subsets<uint32_t>;
template class Subsets<uint64_t>;

template <class bit_t>
SubsetsIterator<bit_t>::SubsetsIterator(int64_t idx) : current_((bit_t)idx) {}

template class SubsetsIterator<uint16_t>;
template class SubsetsIterator<uint32_t>;
template class SubsetsIterator<uint64_t>;

} // namespace xdiag::combinatorics
