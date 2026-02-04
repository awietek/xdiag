// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::combinatorics {

template <typename bit_t> class SubsetsIterator;

// Subsets
template <typename bit_tt> class Subsets {
public:
  using bit_t = bit_tt;
  using iterator_t = SubsetsIterator<bit_t>;

  Subsets() = default;
  explicit Subsets(int64_t n);

  int64_t n() const;
  int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;

private:
  int64_t n_;
  int64_t size_;
  iterator_t begin_, end_;
};

// SubsetsIterator
template <typename bit_t> class SubsetsIterator {
public:
  SubsetsIterator() = default;
  SubsetsIterator(int64_t idx);

  inline bool operator!=(const SubsetsIterator<bit_t> &rhs) const {
    return current_ != rhs.current_;
  }
  inline SubsetsIterator &operator++() {
    ++current_;
    return *this;
  }
  inline bit_t operator*() const { return current_; }

private:
  bit_t current_;
};

} // namespace xdiag::combinatorics
