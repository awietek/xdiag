#pragma once

#include <vector>

#include <xdiag/combinatorics/combinations.h>
#include <xdiag/combinatorics/lin_table.h>
#include <xdiag/combinatorics/subsets.h>
#include <xdiag/common.h>
#include <xdiag/symmetries/permutation_group.h>

namespace xdiag::combinatorics {

template <typename bit_t = std_bit_t> class FermiTableSubsets {
public:
  FermiTableSubsets() = default;
  FermiTableSubsets(int64_t n_sites, PermutationGroup const &group);
  inline bool sign(int64_t sym, bit_t state) const {
    return table_[(sym << n_sites_) | (int64_t)state];
  }
  bool operator==(FermiTableSubsets const &rhs) const;
  bool operator!=(FermiTableSubsets const &rhs) const;

private:
  int64_t n_sites_;
  std::vector<bool> table_;
};

template <typename bit_t = std_bit_t> class FermiTableCombinations {
public:
  FermiTableCombinations() = default;
  FermiTableCombinations(int64_t n_sites, int64_t n_par,
                         PermutationGroup const &group);
  inline bool sign(int64_t sym, bit_t state) const {
    return table_[sym * raw_size_ + lin_table_.index(state)];
  }
  bool operator==(FermiTableCombinations const &rhs) const;
  bool operator!=(FermiTableCombinations const &rhs) const;

private:
  int64_t raw_size_;
  combinatorics::LinTable<bit_t> lin_table_;
  std::vector<bool> table_;
};

} // namespace xdiag::combinatorics
