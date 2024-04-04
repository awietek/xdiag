#pragma once

#include <vector>

#include <xdiag/bits/bitops.h>
#include <xdiag/combinatorics/combinations.h>
#include <xdiag/combinatorics/combinations_index.h>
#include <xdiag/common.h>

namespace xdiag::combinatorics {

template <class bit_t> class LinTable {
public:
  LinTable() = default;
  LinTable(int64_t n, int64_t k);

  inline int64_t n() const { return n_; }
  inline int64_t k() const { return k_; }

  inline int64_t index(bit_t bits) const {
    return left_indices_[bits >> n_right_] +
           right_indices_[bits::gbits(bits, 0, n_right_)];
  }

  inline combinatorics::Combinations<bit_t> states() const {
    return combinatorics::Combinations<bit_t>(n_, k_);
  }

  inline combinatorics::CombinationsIndex<bit_t> states_indices() const {
    return combinatorics::CombinationsIndex<bit_t>(n_, k_);
  }
  inline int64_t size() const { return size_; }

  bool operator==(LinTable<bit_t> const &rhs) const;
  bool operator!=(LinTable<bit_t> const &rhs) const;

private:
  int64_t n_, k_;
  int64_t n_left_;
  int64_t n_right_;

  int64_t left_table_size_;
  int64_t right_table_size_;

  std::vector<int64_t> left_indices_;
  std::vector<int64_t> right_indices_;

  int64_t size_;
};

} // namespace xdiag::combinatorics
