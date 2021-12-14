#pragma once

#include <hydra/common.h>

#include <hydra/combinatorics/bit_patterns.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/combinations_index.h>

#include <hydra/indexing/lintable.h>

namespace hydra::indexing {

template <class bit_t> class CombinationsIndexing {
public:
  CombinationsIndexing() = default;
  CombinationsIndexing(int n, int k);

  inline idx_t index(bit_t state) const { return lin_table_.index(state); }
  inline bit_t operator[](idx_t idx) const {
    return combinatorics::get_nth_pattern<bit_t>(idx, n_, k_);
  }
  inline idx_t size() const { return size_; }

  inline combinatorics::Combinations<bit_t> states() const {
    return combinatorics::Combinations<bit_t>(n_, k_);
  }
  inline combinatorics::CombinationsIndex<bit_t> states_indices() const {
    return combinatorics::CombinationsIndex<bit_t>(n_, k_);
  }

private:
  int n_, k_;
  idx_t size_;
  LinTable<bit_t> lin_table_;
};

} // namespace hydra::indexing
