#pragma once

#include <hydra/common.h>

namespace hydra {

template <class bit_t = std_bit_t> class SubsetsIterator;

// Subsets
template <class bit_t = std_bit_t> class Subsets {
public:
  using iterator_t = SubsetsIterator<bit_t>;

  Subsets() = default;
  explicit Subsets(int n);

  int n() const { return n_; }
  idx_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int n_, k_;
  idx_t size_;
  iterator_t begin_, end_;
};

// SubsetsIterator
template <class bit_t> class SubsetsIterator {
public:

  SubsetsIterator() = default;
  SubsetsIterator(bit_t state);

  inline bool operator==(const SubsetsIterator<bit_t> &rhs) const {
    return current_ == rhs.current_;
  }
  inline bool operator!=(const SubsetsIterator<bit_t> &rhs) const {
    return !operator==(rhs);
  }
  inline SubsetsIterator &operator++() {
    ++current_;
    return *this;
  }
  inline bit_t operator*() const { return current_; }

private:
  bit_t current_;
};

} // namespace hydra
