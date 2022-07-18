#include "permutation_group_lookup.h"

#include <limits>

#include <hydra/symmetries/fermi_sign.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/symmetry_operations.h>

namespace hydra {

template <typename bit_t, int n_sublattice>
GroupActionSublattice<bit_t>::GroupActionSublattice(
    PermutationGroup const &permutation_group)
    : n_sites_(permutation_group.n_sites()),
      n_symmetries_(permutation_group.n_symmetries()),
      permutation_group_(permutation_group), len_word_(n_sites_ / n_sublattice),
      wordmask_((bit_t(1) << len_word_) - 1),
      size_sublat_tables_((idx_t)1 << len_word_) {

  for (int n = 0; n < n_sublat; ++n)
    RepAndSyms.resize(size_sublat_tables_);

  // Loop over all partial words
  for (half_bit_t subword = 0; subword < (half_bit_t)size_sublat_tables_;
       ++subword) {
    
    std::pair<half_bit_t, gsl::span<int const>> RepAndSym[n_sublat];
    for (int n = 0; n < n_sublat; ++n) {
      RepAndSym.first = std::numeric_limits<half_bit_t>::max();
    }

    

    
  }
}

template <typename bit_t>
bool GroupActionSublattice<bit_t>::operator==(
    GroupActionSublattice const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_symmetries_ == rhs.n_symmetries_) &&
         (permutation_group_ == rhs.permutation_group_);
}
template <typename bit_t>
bool GroupActionSublattice<bit_t>::operator!=(
    GroupActionSublattice const &rhs) const {
  return !operator==(rhs);
}

template class GroupActionSublattice<uint16_t>;
template class GroupActionSublattice<uint32_t>;
template class GroupActionSublattice<uint64_t>;

} // namespace hydra
