#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/subsets.h>

#include <hydra/indexing/lintable.h>

namespace hydra::indexing {
template <typename bit_t> class SpinhalfIndexingNoSz {
public:
  SpinhalfIndexingNoSz() = default;
  SpinhalfIndexingNoSz(int n_sites);

  inline int n_sites() const { return n_sites_; }
  inline idx_t size() const { return size_; }
  inline idx_t index(bit_t spins) const { return (idx_t)spins; }

  inline combinatorics::Subsets<bit_t> states() const {
    return combinatorics::Subsets<bit_t>(n_sites_);
  }
  
private:
  int n_sites_;
  idx_t size_;
};

} // namespace hydra::indexing
