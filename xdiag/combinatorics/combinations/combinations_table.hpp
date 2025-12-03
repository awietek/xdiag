// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>

namespace xdiag::combinatorics {

template <class bit_tt> class CombinationsTable {
public:
  using bit_t = bit_tt;

  CombinationsTable() = default;
  CombinationsTable(int64_t n, int64_t k);

  inline int64_t index(bit_t state) const { return lin_table_.index(state); }
  bit_t operator[](int64_t idx) const;
  int64_t size() const;
  Combinations<bit_t> states() const;

  bool operator==(CombinationsTable<bit_t> const &rhs) const;
  bool operator!=(CombinationsTable<bit_t> const &rhs) const;

private:
  int64_t n_, k_;
  int64_t size_;
  LinTable<bit_t> lin_table_;
};

template <class bit_t>
CombinationsTable<bit_t> table(Combinations<bit_t> const &states);

} // namespace xdiag::combinatorics
