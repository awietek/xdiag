#pragma once

#include <utility>
#include <vector>

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra::indexing {

template <class bit_t, class GroupAction> class IndexingSymmetricFermionic {
public:
  IndexingSymmetricFermionic() = default;
  IndexingSymmetricFermionic(int n_sites, int n_up,
                             PermutationGroup permutation_group,
                             Representation irrep);

  inline idx_t size() const { return size_; }
  inline bit_t state(idx_t idx) const { return states_[idx]; }
  inline idx_t index(bit_t state) const {
    return index_of_raw_state_[lin_table_.index(state)];
  }
  inline complex norm(bit_t state) const {
    return norm_of_raw_state_[lin_table_.index(state)];
  }

private:
  LinTable<bit_t> lin_table_;
  std::vector<bit_t> states_;
  std::vector<idx_t> index_of_raw_state_;
  std::vector<complex> norm_of_raw_state_;

  idx_t size_;
};

} // namespace hydra::indexing
