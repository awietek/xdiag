#pragma once

#include <hydra/common.h>
#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/bit_patterns.h>

namespace hydra {

template <class bit_t = std_bit_t> class CombinationsIterator;

// Combinations
template <class bit_t = std_bit_t> class Combinations {
public:
  using iterator_t = CombinationsIterator<bit_t>;

  Combinations() = default;
  Combinations(int n, int k);

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

// CombinationsIterator
template <class bit_t> class CombinationsIterator {
public:

  CombinationsIterator() = default;
  CombinationsIterator(bit_t state);

  inline bool operator==(CombinationsIterator<bit_t> const &rhs) const {
    return current_ == rhs.current_;
  }
  inline bool operator!=(CombinationsIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }
  inline CombinationsIterator &operator++() {
    current_ = combinatorics::get_next_pattern(current_);
    return *this;
  }
  inline bit_t operator*() const { return current_; }

private:
  bit_t current_;
};

} // namespace hydra
