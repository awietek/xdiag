// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/combinations_index.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/common.hpp>

namespace xdiag::combinatorics {

template <class bit_t> class CombinationsIndexing {
public:
  CombinationsIndexing() = default;
  CombinationsIndexing(int64_t n, int64_t k);

  inline int64_t index(bit_t state) const { return lin_table_.index(state); }
  inline bit_t operator[](int64_t idx) const {
    return combinatorics::get_nth_pattern<bit_t>(idx, n_, k_);
  }
  inline int64_t size() const { return size_; }

  inline combinatorics::Combinations<bit_t> states() const {
    return combinatorics::Combinations<bit_t>(n_, k_);
  }
  inline combinatorics::CombinationsIndex<bit_t> states_indices() const {
    return combinatorics::CombinationsIndex<bit_t>(n_, k_);
  }

private:
  int64_t n_, k_;
  int64_t size_;
  LinTable<bit_t> lin_table_;
};

} // namespace xdiag::combinatorics
