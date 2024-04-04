#include "lin_table.hpp"

#include <xdiag/bits/bitops.hpp>
#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations.hpp>

namespace xdiag::combinatorics {

template <class bit_t>
LinTable<bit_t>::LinTable(int64_t n, int64_t k)
    : n_(n), k_(k), n_left_(ceil(n / 2.0)), n_right_(n - n_left_),
      left_table_size_((int64_t)1 << n_left_),
      right_table_size_((int64_t)1 << n_right_),
      left_indices_(left_table_size_, 0), right_indices_(right_table_size_, 0),
      size_(combinatorics::binomial(n, k)) {
  using bits::popcnt;

  // Fill offsets on left indices
  for (bit_t left = 0; left < (bit_t)left_table_size_; ++left) {

    if (left == 0) {
      left_indices_[left] = 0;
    } else {
      left_indices_[left] = left_indices_[left - 1] +
                            binomial(n_right_, k - popcnt(bit_t(left - 1)));
    }
  }

  // Fill indices for right combinations
  for (int64_t k_right = 0; k_right <= n_right_; ++k_right) {
    int64_t idx = 0;
    for (bit_t bits : Combinations<bit_t>(n_right_, k_right))
      right_indices_[bits] = idx++;
  }
}

template <class bit_t>
bool LinTable<bit_t>::operator==(LinTable<bit_t> const &rhs) const {
  return (n_ == rhs.n_) && (k_ == rhs.k_);
}

template <class bit_t>
bool LinTable<bit_t>::operator!=(LinTable<bit_t> const &rhs) const {
  return !operator==(rhs);
}

template class LinTable<uint16_t>;
template class LinTable<uint32_t>;
template class LinTable<uint64_t>;

} // namespace xdiag::combinatorics
