#pragma once

#include <hydra/common.h>
#include <hydra/combinatorics/subsets.h>

namespace hydra::combinatorics {

template <class bit_t> class SubsetsIndexIterator;

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
  SubsetsIndexIterator(idx_t idx);

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

#ifdef HYDRA_ENABLE_OPENMP
// SubsetsIndexThread
template <typename bit_t = std_bit_t> class SubsetsIndexThread {
public:
  using iterator_t = SubsetsIndexIterator<bit_t>;

  SubsetsIndexThread() = default;
  SubsetsIndexThread(int n);

  int n() const { return n_; }
  idx_t size() const { return size_; };
  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

private:
  int n_;
  idx_t size_;
  iterator_t begin_, end_;
};

template <typename bit_t>
inline SubsetsIndexThread<bit_t> ThreadStatesIndex(Subsets<bit_t> const &si) {
  return SubsetsIndexThread(si.n());
}

template <typename bit_t>
inline SubsetsIndexThread<bit_t>
ThreadStatesIndex(SubsetsIndex<bit_t> const &si) {
  return SubsetsIndexThread(si.n());
}
#endif

} // namespace hydra::combinatorics
