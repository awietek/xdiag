#pragma once

#include <vector>

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/combinatorics/subsets.h>
#include <hydra/combinatorics/subsets_index.h>

namespace hydra::indexing {

template <class bit_t = std_bit_t> class LinTableRaw {
public:
  LinTableRaw() = default;
  LinTableRaw(int n);

  inline idx_t index(bit_t bits) const { return (idx_t)bits; }

  inline combinatorics::Subsets<bit_t> states() const {
    return combinatorics::Subsets<bit_t>(n_);
  }

  inline combinatorics::SubsetsIndex<bit_t> states_indices() const {
    return combinatorics::SubsetsIndex<bit_t>(n_);
  }
  inline idx_t size() const { return size_; }


private:
  int n_;
  idx_t size_;
};

// Comment: Rather dumb class, drop in replacement for real LinTable without
// conservation

} // namespace hydra::indexing
