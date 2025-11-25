// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/common.hpp>

namespace xdiag::combinatorics {

template <class bit_t> class SubsetsIndexIterator;

// SubsetsIndex
template <class bit_type> class SubsetsIndex {
public:
  using bit_t = bit_type;
  using iterator_t = SubsetsIndexIterator<bit_t>;

  SubsetsIndex() = default;
  explicit SubsetsIndex(int64_t n);

  int64_t n() const { return n_; }
  int64_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int64_t n_, k_;
  int64_t size_;
  iterator_t begin_, end_;
};

// SubsetsIndexIterator
template <class bit_t> class SubsetsIndexIterator {
public:
  SubsetsIndexIterator() = default;
  SubsetsIndexIterator(int64_t idx);

  inline bool operator==(const SubsetsIndexIterator<bit_t> &rhs) const {

    return current_ == rhs.current_;
  }
  inline bool operator!=(const SubsetsIndexIterator<bit_t> &rhs) const {
    return !operator==(rhs);
  }
  inline SubsetsIndexIterator &operator++() {
    ++current_;
    return *this;
  }
  inline std::pair<bit_t, int64_t> operator*() const {
    return {current_, (int64_t)current_};
  }

private:
  bit_t current_;
};

#ifdef _OPENMP
// SubsetsIndexThread
template <typename bit_t> class SubsetsIndexThread {
public:
  using iterator_t = SubsetsIndexIterator<bit_t>;

  SubsetsIndexThread() = default;
  SubsetsIndexThread(int64_t n);

  int64_t n() const { return n_; }
  int64_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int64_t n_;
  int64_t size_;
  iterator_t begin_, end_;
};

template <typename bit_t>
inline SubsetsIndexThread<bit_t> ThreadStatesIndex(Subsets<bit_t> const &si) {
  return SubsetsIndexThread<bit_t>(si.n());
}

template <typename bit_t>
inline SubsetsIndexThread<bit_t>
ThreadStatesIndex(SubsetsIndex<bit_t> const &si) {
  return SubsetsIndexThread<bit_t>(si.n());
}
#endif

} // namespace xdiag::combinatorics
