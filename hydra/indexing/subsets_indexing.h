#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/bit_patterns.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/combinatorics/subsets_index.h>

namespace hydra::indexing {

template <class bit_t> class SubsetsIndexing {
public:
  SubsetsIndexing() = default;
  SubsetsIndexing(int n);

  inline idx_t index(bit_t state) const { return (idx_t)state; }
  inline bit_t operator[](idx_t idx) const { return (bit_t)idx; }
  inline idx_t size() const { return size_; }

  inline combinatorics::Subsets<bit_t> states() const {
    return combinatorics::Subsets<bit_t>(n_);
  }
  inline combinatorics::SubsetsIndex<bit_t> states_indices() const {
    return combinatorics::SubsetsIndex<bit_t>(n_);
  }

private:
  int n_;
  idx_t size_;
};

} // namespace hydra::indexing
