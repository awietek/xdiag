// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/combinatorics/subsets_index.hpp>
#include <xdiag/common.hpp>

namespace xdiag::combinatorics {

template <class bit_t> class SubsetsIndexing {
public:
  SubsetsIndexing() = default;
  SubsetsIndexing(int64_t n);

  inline int64_t index(bit_t state) const { return (int64_t)state; }
  inline bit_t operator[](int64_t idx) const { return (bit_t)idx; }
  inline int64_t size() const { return size_; }

  inline combinatorics::Subsets<bit_t> states() const {
    return combinatorics::Subsets<bit_t>(n_);
  }
  inline combinatorics::SubsetsIndex<bit_t> states_indices() const {
    return combinatorics::SubsetsIndex<bit_t>(n_);
  }

private:
  int64_t n_;
  int64_t size_;
};

} // namespace xdiag::combinatorics
