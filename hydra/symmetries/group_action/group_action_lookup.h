#pragma once
#include <tuple>

#include <hydra/common.h>
#include <hydra/symmetries/permutation_group.h>

namespace hydra {

template <class bit_t> class GroupActionLookup {
public:
  GroupActionLookup() = default;
  explicit GroupActionLookup(PermutationGroup const &permutation_group);

  inline int n_sites() const { return n_sites_; }
  inline int n_symmetries() const { return n_symmetries_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }

  inline bit_t apply(int sym, bit_t state) const {
    // return table_prefix_[sym * prefix_size_ + (state >> n_postfix_bits_)] |
    //        table_postfix_[sym * postfix_size_ + (state & postfix_mask_)];
    return table_prefix_[(sym << n_prefix_bits_) | (state >> n_postfix_bits_)] |
           table_postfix_[(sym << n_postfix_bits_) | (state & postfix_mask_)];
  }

  bool operator==(GroupActionLookup const &rhs) const;
  bool operator!=(GroupActionLookup const &rhs) const;

private:
  int n_sites_;
  int n_symmetries_;
  PermutationGroup permutation_group_;

  int n_prefix_bits_;
  int n_postfix_bits_;
  bit_t postfix_mask_;

  idx_t prefix_size_;
  idx_t postfix_size_;

  std::vector<bit_t> table_prefix_;
  std::vector<bit_t> table_postfix_;
};

} // namespace hydra
