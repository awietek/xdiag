#include "group_action_lookup.hpp"

#include <xdiag/symmetries/group_action/group_action.hpp>
#include <xdiag/symmetries/operations/fermi_sign.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>

namespace xdiag {

template <typename bit_t>
GroupActionLookup<bit_t>::GroupActionLookup(
    PermutationGroup const &permutation_group)
    : n_sites_(permutation_group.n_sites()),
      n_symmetries_(permutation_group.n_symmetries()),
      permutation_group_(permutation_group), n_prefix_bits_(n_sites_ / 2),
      n_postfix_bits_(n_sites_ - n_prefix_bits_),
      postfix_mask_(((bit_t)1 << n_postfix_bits_) - 1),
      prefix_size_(pow(2, n_prefix_bits_)),
      postfix_size_(pow(2, n_postfix_bits_)),
      table_prefix_(n_symmetries_ * prefix_size_, 0),
      table_postfix_(n_symmetries_ * postfix_size_, 0) {

  auto action = GroupAction(permutation_group);

  int64_t idx = 0;
  // Fill prefix table with translated prefix states
  for (int sym = 0; sym < n_symmetries_; ++sym) {
    for (bit_t state = 0; state < (bit_t)prefix_size_; ++state) {
      table_prefix_[idx++] =
          action.apply(sym, (bit_t)(state << n_postfix_bits_));
    }
  }
  assert(idx == (int64_t)table_prefix_.size());

  // Fill postfix table
  idx = 0;
  for (int sym = 0; sym < n_symmetries_; ++sym) {
    for (bit_t state = 0; state < (bit_t)postfix_size_; ++state) {
      table_postfix_[idx++] = action.apply(sym, state);
    }
  }

  assert(idx == (int64_t)table_postfix_.size());
}

template <typename bit_t>
bool GroupActionLookup<bit_t>::operator==(GroupActionLookup const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_symmetries_ == rhs.n_symmetries_) &&
         (permutation_group_ == rhs.permutation_group_);
}
template <typename bit_t>
bool GroupActionLookup<bit_t>::operator!=(GroupActionLookup const &rhs) const {
  return !operator==(rhs);
}

template class GroupActionLookup<uint16_t>;
template class GroupActionLookup<uint32_t>;
template class GroupActionLookup<uint64_t>;

} // namespace xdiag
