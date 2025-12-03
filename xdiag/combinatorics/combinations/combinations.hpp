// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/bit_patterns.hpp>

namespace xdiag::combinatorics {

template <typename bit_t> class CombinationsIterator;

template <typename bit_tt> class Combinations {
public:
  using bit_t = bit_tt;
  using iterator_t = CombinationsIterator<bit_t>;

  Combinations() = default;
  Combinations(int64_t n, int64_t k);

  int64_t n() const;
  int64_t k() const;

  int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(Combinations<bit_t> const &rhs) const;
  bool operator!=(Combinations<bit_t> const &rhs) const;

private:
  int64_t n_, k_;
  int64_t size_;
  iterator_t begin_, end_;
};

template <typename bit_t> class CombinationsIterator {
public:
  CombinationsIterator() = default;
  CombinationsIterator(int64_t n, int64_t k, int64_t idx);
  inline bool operator!=(CombinationsIterator<bit_t> const &rhs) const {
    return idx_ != rhs.idx_;
  }
  inline CombinationsIterator &operator++() {
    current_ = combinatorics::get_next_pattern(current_);
    ++idx_;
    return *this;
  }
  inline bit_t operator*() const { return current_; }

private:
  bit_t current_;
  int64_t idx_;
};

} // namespace xdiag::combinatorics
