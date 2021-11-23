#pragma once

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>

namespace hydra::indexing {
template <typename bit_t> class SpinhalfIndexing {
public:
  SpinhalfIndexing() = default;
  SpinhalfIndexing(int n_sites, int nup);

  inline int n_sites() const { return n_sites_; }
  inline int n_up() const { return n_up_; }
  inline idx_t size() const { return size_; }
  inline idx_t index(bit_t spins) const { return lintable_.index(spins); }

private:
  int n_sites_;
  int n_up_;
  indexing::LinTable<bit_t> lintable_;
  idx_t size_;
};

} // namespace hydra::indexing
