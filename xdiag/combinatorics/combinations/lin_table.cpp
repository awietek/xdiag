// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "lin_table.hpp"

#include <cmath>

#include <xdiag/bits/popcount.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/math/binomial.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
LinTable<bit_t>::LinTable(int64_t n, int64_t k) try
    : combinations_(n, k), n_left_(ceil(n / 2.0)), n_right_(n - n_left_),
      left_table_size_((int64_t)1 << n_left_),
      right_table_size_((int64_t)1 << n_right_),
      left_indices_(left_table_size_, 0), right_indices_(right_table_size_, 0) {

  // Fill offsets on left indices
  for (bit_t left = 0; left < (bit_t)left_table_size_; ++left) {

    if (left == 0) {
      left_indices_[left] = 0;
    } else {
      left_indices_[left] =
          left_indices_[left - 1] +
          math::binomial(n_right_, k - bits::popcount(bit_t(left - 1)));
    }
  }

  // Fill indices for right combinations
  for (int64_t k_right = 0; k_right <= n_right_; ++k_right) {
    int64_t idx = 0;
    for (bit_t bits : Combinations<bit_t>(n_right_, k_right)) {
      right_indices_[bits] = idx++;
    }
  }
}
XDIAG_CATCH

template <class bit_t> int64_t LinTable<bit_t>::n() const {
  return combinations_.n();
}

template <class bit_t> int64_t LinTable<bit_t>::k() const {
  return combinations_.k();
}

template <class bit_t> int64_t LinTable<bit_t>::size() const {
  return combinations_.size();
}

template <class bit_t> int64_t LinTable<bit_t>::bitwidth() const {
  return combinations_.n();
}

template <class bit_t> bit_t LinTable<bit_t>::operator[](int64_t idx) const {
  return combinations_[idx];
}

template <class bit_t>
CombinationsIterator<bit_t> LinTable<bit_t>::begin() const {
  return CombinationsIterator<bit_t>(combinations_.n(), combinations_.k(), 0);
}

template <class bit_t>
CombinationsIterator<bit_t> LinTable<bit_t>::end() const {
  return CombinationsIterator<bit_t>(combinations_.n(), combinations_.k(),
                                     combinations_.size());
}

template <class bit_t>
bool LinTable<bit_t>::operator==(LinTable<bit_t> const &rhs) const {
  return combinations_ == rhs.combinations_;
}

template <class bit_t>
bool LinTable<bit_t>::operator!=(LinTable<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class LinTable<uint32_t>;
template class LinTable<uint64_t>;

} // namespace xdiag::combinatorics
