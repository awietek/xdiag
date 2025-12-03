// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "combinations.hpp"

#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
Combinations<bit_t>::Combinations(int64_t n, int64_t k) try
    : n_(n), k_(k), size_(combinatorics::binomial(n, k)) {
  if (k > n) {
    XDIAG_THROW("Error constructing Combinations: k>n");
  } else if (k < 0) {
    XDIAG_THROW("Error constructing Combinations: k<0");
  } else if (n < 0) {
    XDIAG_THROW("Error constructing Combinations: n<0");
  } else {
    begin_ = CombinationsIterator<bit_t>(n, k, (int64_t)0);
    end_ = CombinationsIterator<bit_t>(n, k, (int64_t)binomial(n, k));
  }
}
XDIAG_CATCH

template <class bit_t> int64_t Combinations<bit_t>::n() const { return n_; }
template <class bit_t> int64_t Combinations<bit_t>::k() const { return k_; }
template <class bit_t> int64_t Combinations<bit_t>::size() const {
  return size_;
};

template <class bit_t>
CombinationsIterator<bit_t> Combinations<bit_t>::begin() const {
  return begin_;
}

template <class bit_t>
CombinationsIterator<bit_t> Combinations<bit_t>::end() const {
  return end_;
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

template <class bit_t>
CombinationsIterator<bit_t>::CombinationsIterator(int64_t n, int64_t k,
                                                  int64_t idx)
    : current_(get_nth_pattern<bit_t>(idx, n, k)), idx_(idx) {}

template class CombinationsIterator<uint16_t>;
template class CombinationsIterator<uint32_t>;
template class CombinationsIterator<uint64_t>;

} // namespace xdiag::combinatorics
