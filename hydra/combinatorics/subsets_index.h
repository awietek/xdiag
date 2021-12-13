#pragma once

#include <hydra/common.h>

namespace hydra::combinatorics {

template <class bit_t = std_bit_t> class SubsetsIndexIterator;

// SubsetsIndex
template <class bit_t = std_bit_t> class SubsetsIndex {
public:
  using iterator_t = SubsetsIndexIterator<bit_t>;

  SubsetsIndex() = default;
  explicit SubsetsIndex(int n);

  int n() const { return n_; }
  idx_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int n_, k_;
  idx_t size_;
  iterator_t begin_, end_;
};

// SubsetsIndexIterator
template <class bit_t> class SubsetsIndexIterator {
public:
  SubsetsIndexIterator() = default;
  SubsetsIndexIterator(bit_t state);

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
  inline std::pair<bit_t, idx_t> operator*() const {
    return {current_, (idx_t)current_};
  }

private:
  bit_t current_;
};

} // namespace hydra::combinatorics
