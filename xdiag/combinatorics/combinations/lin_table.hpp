// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>
#include <xdiag/bits/bitops.hpp>

namespace xdiag::combinatorics {

template <class bit_tt> class LinTable {
public:
  using bit_t = bit_tt;
  LinTable() = default;
  LinTable(int64_t n, int64_t k);
  inline int64_t index(bit_t bits) const {
    return left_indices_[bits >> n_right_] +
           right_indices_[bits::gbits(bits, 0, n_right_)];
  }
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
