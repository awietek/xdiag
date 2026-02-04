// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/combinatorics/subsets/subsets.hpp>

namespace xdiag::combinatorics {

template <typename bit_tt> class SubsetsTable {
public:
  using bit_t = bit_tt;

  SubsetsTable() = default;
  SubsetsTable(int64_t n);

  inline int64_t index(bit_t state) const { return (int64_t)state; }
  inline bit_t operator[](int64_t idx) const { return (bit_t)idx; }
  int64_t size() const;
  Subsets<bit_t> states() const;

private:
  int64_t n_;
  int64_t size_;
};

} // namespace xdiag::combinatorics
