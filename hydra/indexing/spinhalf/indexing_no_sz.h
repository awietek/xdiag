#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/subsets_index.h>

namespace hydra::indexing::spinhalf {

template <typename bit_t> class IndexingNoSz {
public:
  using iterator_t = combinatorics::SubsetsIndexIterator<bit_t>;

  IndexingNoSz() = default;
  IndexingNoSz(int n_sites);

  inline int n_sites() const { return n_sites_; }
  inline idx_t size() const { return size_; }
  inline idx_t index(bit_t spins) const { return (idx_t)spins; }

  iterator_t begin() const { return begin_; }
  iterator_t end() const { return end_; }

#ifdef HYDRA_ENABLE_OPENMP
  iterator_t thread_begin() const;
  iterator_t thread_end() const;
#endif

private:
  int n_sites_;
  idx_t size_;
  iterator_t begin_, end_;
};

} // namespace hydra::indexing::spinhalf
