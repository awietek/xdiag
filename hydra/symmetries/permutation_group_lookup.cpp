#include "permutation_group_lookup.h"

#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/symmetry_utils.h>

namespace hydra {

template <class bit_t>
PermutationGroupLookup<bit_t>::PermutationGroupLookup(
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

  auto action = PermutationGroupAction(permutation_group);

  idx_t idx = 0;
  // Fill prefix table with translated prefix states
  for (int sym = 0; sym < n_symmetries_; ++sym) {
    for (bit_t state = 0; state < prefix_size_; ++state) {
      table_prefix_[idx++] =
          action.apply(sym, (bit_t)(state << n_postfix_bits_));
    }
  }
  assert(idx == table_prefix_.size());

  // Fill postfix table
  idx = 0;
  for (int sym = 0; sym < n_symmetries_; ++sym) {
    for (bit_t state = 0; state < postfix_size_; ++state) {
      table_postfix_[idx++] = action.apply(sym, state);
    }
  }

  assert(idx == table_postfix_.size());
}

template <class bit_t>
bool PermutationGroupLookup<bit_t>::operator==(
    PermutationGroupLookup const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_symmetries_ == rhs.n_symmetries_) &&
         (permutation_group_ == rhs.permutation_group_);
}
template <class bit_t>
bool PermutationGroupLookup<bit_t>::operator!=(
    PermutationGroupLookup const &rhs) const {
  return !operator==(rhs);
}

  template class PermutationGroupLookup<uint16_t>;
template class PermutationGroupLookup<uint32_t>;
template class PermutationGroupLookup<uint64_t>;

} // namespace hydra
