#pragma once

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/common.hpp>

namespace xdiag::combinatorics {

template <typename bit_t> class CombinationsIterator;

// Combinations
template <typename bit_t> class Combinations {
public:
  using iterator_t = CombinationsIterator<bit_t>;

  Combinations() = default;
  Combinations(int64_t n, int64_t k);

  inline int64_t n() const { return n_; }
  inline int64_t k() const { return k_; }
  inline int64_t size() const { return size_; };
  inline iterator_t begin() const { return begin_; }
  inline iterator_t end() const { return end_; }

private:
  int64_t n_, k_;
  int64_t size_;
  iterator_t begin_, end_;
};

// CombinationsIterator
template <typename bit_t> class CombinationsIterator {
public:
  CombinationsIterator() = default;
  CombinationsIterator(int64_t n, int64_t k, int64_t idx);

  inline bool operator==(CombinationsIterator<bit_t> const &rhs) const {
    return idx_ == rhs.idx_;
  }
  inline bool operator!=(CombinationsIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
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

#ifdef _OPENMP
// CombinationsThread
template <typename bit_t> class CombinationsThread {
public:
  using iterator_t = CombinationsIterator<bit_t>;

  CombinationsThread() = default;
  CombinationsThread(int64_t n, int64_t k);

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
inline CombinationsThread<bit_t> ThreadStates(Combinations<bit_t> const &si) {
  return CombinationsThread<bit_t>(si.n(), si.k());
}

#endif

} // namespace xdiag::combinatorics
