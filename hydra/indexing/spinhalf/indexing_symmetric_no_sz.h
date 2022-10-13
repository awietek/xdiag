#pragma once

#include <extern/gsl/span>
#include <utility>
#include <vector>

#include <hydra/common.h>
#include <hydra/indexing/subsets_indexing.h>
#include <hydra/indexing/spinhalf/symmetric_iterator.h>
#include <hydra/symmetries/group_action/group_action_lookup.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra::indexing::spinhalf {

template <typename bit_t> class IndexingSymmetricNoSz {
public:
  using iterator_t = SymmetricIterator<bit_t>;
  using span_size_t = gsl::span<int const>::size_type;

  IndexingSymmetricNoSz() = default;
  IndexingSymmetricNoSz(int n_sites, PermutationGroup permutation_group,
                        Representation irrep);
  inline int n_sites() const { return n_sites_; }

  GroupActionLookup<bit_t> const &group_action() const { return group_action_; }
  Representation const &irrep() const { return irrep_; }
  inline idx_t size() const { return size_; }

  inline bit_t state(idx_t idx) const { return reps_[idx]; }
  inline double norm(idx_t idx) const { return norms_[idx]; }

  inline idx_t index(bit_t state) const {
    return index_for_rep_[subsets_indexing_.index(state)];
  }
  inline bit_t representative(bit_t raw_state) const {
    return reps_[index(raw_state)];
  }

  inline std::pair<idx_t, int> index_sym(bit_t raw_state) const {
    idx_t raw_idx = subsets_indexing_.index(raw_state);
    idx_t index = index_for_rep_[raw_idx];
    if (index == invalid_index) {
      return {invalid_index, 0};
    }
    idx_t start = sym_limits_for_rep_[raw_idx].first;
    return {index, syms_[start]};
  }

  inline std::pair<idx_t, gsl::span<int const>>
  index_syms(bit_t raw_state) const {
    idx_t raw_idx = subsets_indexing_.index(raw_state);
    idx_t index = index_for_rep_[raw_idx];
    auto [start, length] = sym_limits_for_rep_[raw_idx];
    return {index, {syms_.data() + start, length}};
  }

  inline iterator_t begin() const { return begin_; }
  inline iterator_t end() const { return end_; }
  
  bool operator==(IndexingSymmetricNoSz const &rhs) const;
  bool operator!=(IndexingSymmetricNoSz const &rhs) const;
  
private:
  int n_sites_;
  GroupActionLookup<bit_t> group_action_;
  Representation irrep_;
  SubsetsIndexing<bit_t> subsets_indexing_;

  std::vector<bit_t> reps_;
  std::vector<idx_t> index_for_rep_;
  std::vector<int> syms_;
  std::vector<std::pair<span_size_t, span_size_t>> sym_limits_for_rep_;
  std::vector<double> norms_;

  idx_t size_;
  iterator_t begin_, end_; 
};

} // namespace hydra::indexing::spinhalf
