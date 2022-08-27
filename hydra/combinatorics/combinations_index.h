#pragma once

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/bit_patterns.h>
#include <hydra/common.h>

namespace hydra::combinatorics {

template <typename bit_t = std_bit_t> class CombinationsIndexIterator;

// CombinationsIndex
template <typename bit_t = std_bit_t> class CombinationsIndex {
public:
  using iterator_t = CombinationsIndexIterator<bit_t>;

  CombinationsIndex() = default;
  CombinationsIndex(int n, int k);

  int n() const { return n_; }
  int k() const { return k_; }
  idx_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int n_, k_;
  idx_t size_;
  iterator_t begin_, end_;
};

// CombinationsIndexIterator
template <typename bit_t> class CombinationsIndexIterator {
public:
  CombinationsIndexIterator() = default;
  CombinationsIndexIterator(int n, int k, idx_t idx);

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
  inline std::pair<bit_t, idx_t> operator*() const { return {current_, idx_}; }
  inline idx_t idx() const { return idx_; }

private:
  bit_t current_;
  idx_t idx_;
};

#ifdef HYDRA_ENABLE_OPENMP
// CombinationsIndexThread
template <typename bit_t = std_bit_t> class CombinationsIndexThread {
public:
  using iterator_t = CombinationsIndexIterator<bit_t>;

  CombinationsIndexThread() = default;
  CombinationsIndexThread(int n, int k);

  int n() const { return n_; }
  int k() const { return k_; }
  idx_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int n_, k_;
  idx_t size_;
  iterator_t begin_, end_;
};

template <typename bit_t>
inline CombinationsIndexThread<bit_t>
ThreadStatesIndex(Combinations<bit_t> const &si) {
  return CombinationsIndexThread(si.n(), si.k());
}

template <typename bit_t>
inline CombinationsIndexThread<bit_t>
ThreadStatesIndex(CombinationsIndex<bit_t> const &si) {
  return CombinationsIndexThread(si.n(), si.k());
}
#endif
} // namespace hydra::combinatorics
