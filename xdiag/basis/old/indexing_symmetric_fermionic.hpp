#pragma once

#include <utility>
#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/indexing/lintable.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::indexing {

template <typename bit_t, class GroupAction> class IndexingSymmetricFermionic {
public:
  IndexingSymmetricFermionic() = default;
  IndexingSymmetricFermionic(int n_sites, int n_fermions,
                             PermutationGroup permutation_group,
                             Representation irrep);

  inline int64_t size() const { return size_; }
  inline bit_t state(int64_t idx) const { return states_[idx]; }
  inline int64_t index(bit_t state) const {
    return index_of_raw_state_[lin_table_.index(state)];
  }
  inline complex norm(bit_t state) const {
    return norm_of_raw_state_[lin_table_.index(state)];
  }

  inline auto begin() const { return states_.begin(); }
  inline auto end() const { return states_.end(); }

private:
  LinTable<bit_t> lin_table_;
  std::vector<bit_t> states_;
  std::vector<int64_t> index_of_raw_state_;
  std::vector<complex> norm_of_raw_state_;

  int64_t size_;
};

} // namespace xdiag::indexing
