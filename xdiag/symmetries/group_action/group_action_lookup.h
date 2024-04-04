#pragma once
#include <tuple>

#include <xdiag/common.h>
#include <xdiag/symmetries/permutation_group.h>

namespace xdiag {

template <class bit_t> class GroupActionLookup {
public:
  GroupActionLookup() = default;
  explicit GroupActionLookup(PermutationGroup const &permutation_group);

  inline int64_t n_sites() const { return n_sites_; }
  inline int64_t n_symmetries() const { return n_symmetries_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }

  inline bit_t apply(int64_t sym, bit_t state) const {
    // return table_prefix_[sym * prefix_size_ + (state >> n_postfix_bits_)] |
    //        table_postfix_[sym * postfix_size_ + (state & postfix_mask_)];
    return table_prefix_[(sym << n_prefix_bits_) | (state >> n_postfix_bits_)] |
      table_postfix_[(sym << n_postfix_bits_) | (state & postfix_mask_)];
  }

  bool operator==(GroupActionLookup const &rhs) const;
  bool operator!=(GroupActionLookup const &rhs) const;

private:
  int64_t n_sites_;
  int64_t n_symmetries_;
  PermutationGroup permutation_group_;

  int64_t n_prefix_bits_;
  int64_t n_postfix_bits_;
  bit_t postfix_mask_;

  int64_t prefix_size_;
  int64_t postfix_size_;

  std::vector<bit_t> table_prefix_;
  std::vector<bit_t> table_postfix_;
};

} // namespace xdiag
