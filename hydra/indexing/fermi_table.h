#pragma once

#include <vector>

#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/common.h>
#include <hydra/indexing/lin_table.h>
#include <hydra/symmetries/permutation_group.h>

namespace hydra::indexing {

template <typename bit_t = std_bit_t> class FermiTableSubsets {
public:
  FermiTableSubsets() = default;
  FermiTableSubsets(int n_sites, PermutationGroup const &group);
  inline bool sign(int sym, bit_t state) const {
    return table_[(sym << n_sites_) | (idx_t)state];
  }
  bool operator==(FermiTableSubsets const &rhs) const;
  bool operator!=(FermiTableSubsets const &rhs) const;

private:
  int n_sites_;
  std::vector<bool> table_;
};

template <typename bit_t = std_bit_t> class FermiTableCombinations {
public:
  FermiTableCombinations() = default;
  FermiTableCombinations(int n_sites, int n_par, PermutationGroup const &group);
  inline bool sign(int sym, bit_t state) const {
    return table_[sym * raw_size_ + lin_table_.index(state)];
  }
  bool operator==(FermiTableCombinations const &rhs) const;
  bool operator!=(FermiTableCombinations const &rhs) const;

private:
  idx_t raw_size_;
  LinTable<bit_t> lin_table_;
  std::vector<bool> table_;
};

} // namespace hydra::indexing
