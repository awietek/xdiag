#pragma once

#include <hydra/common.h>

namespace hydra::combinatorics {

template <class bit_t> class SubsetsIterator;

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
  SubsetsIterator(idx_t idx);

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

#ifdef HYDRA_ENABLE_OPENMP
// SubsetsThread
template <class bit_t = std_bit_t> class SubsetsThread {
public:
  using iterator_t = SubsetsIterator<bit_t>;

  SubsetsThread() = default;
  explicit SubsetsThread(int n);

  int n() const { return n_; }
  idx_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int n_, k_;
  idx_t size_;
  iterator_t begin_, end_;
};

template <typename bit_t>
SubsetsThread<bit_t> ThreadStates(Subsets<bit_t> const &si) {
  return SubsetsThread<bit_t>(si.n());
}
#endif

} // namespace hydra::combinatorics
