#pragma once

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/common.hpp>

namespace xdiag::combinatorics {

template <typename bit_t> class CombinationsIndexIterator;

// CombinationsIndex
template <typename bit_t> class CombinationsIndex {
public:
  using iterator_t = CombinationsIndexIterator<bit_t>;

  CombinationsIndex() = default;
  CombinationsIndex(int64_t n, int64_t k);

  int64_t n() const { return n_; }
  int64_t k() const { return k_; }
  int64_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int64_t n_, k_;
  int64_t size_;
  iterator_t begin_, end_;
};

// CombinationsIndexIterator
template <typename bit_t> class CombinationsIndexIterator {
public:
  CombinationsIndexIterator() = default;
  CombinationsIndexIterator(int64_t n, int64_t k, int64_t idx);

  inline bool operator==(CombinationsIndexIterator<bit_t> const &rhs) const {
    return idx_ == rhs.idx_;
  }
  inline bool operator!=(CombinationsIndexIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }
  inline CombinationsIndexIterator &operator++() {
    current_ = combinatorics::get_next_pattern(current_);
    ++idx_;
    return *this;
  }
  inline std::pair<bit_t, int64_t> operator*() const { return {current_, idx_}; }
  inline int64_t idx() const { return idx_; }

private:
  bit_t current_;
  int64_t idx_;
};

#ifdef _OPENMP
// CombinationsIndexThread
template <typename bit_t> class CombinationsIndexThread {
public:
  using iterator_t = CombinationsIndexIterator<bit_t>;

  CombinationsIndexThread() = default;
  CombinationsIndexThread(int64_t n, int64_t k);

  int64_t n() const { return n_; }
  int64_t k() const { return k_; }
  int64_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int64_t n_, k_;
  int64_t size_;
  iterator_t begin_, end_;
};

template <typename bit_t>
inline CombinationsIndexThread<bit_t>
ThreadStatesIndex(Combinations<bit_t> const &si) {
  return CombinationsIndexThread<bit_t>(si.n(), si.k());
}

template <typename bit_t>
inline CombinationsIndexThread<bit_t>
ThreadStatesIndex(CombinationsIndex<bit_t> const &si) {
  return CombinationsIndexThread<bit_t>(si.n(), si.k());
}
#endif
} // namespace xdiag::combinatorics
