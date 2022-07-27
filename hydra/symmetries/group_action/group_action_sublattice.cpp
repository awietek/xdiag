#include "group_action_sublattice.h"

#include <limits>

#include <hydra/symmetries/group_action/group_action.h>
#include <hydra/symmetries/group_action/sublattice_stability.h>

#include <hydra/utils/logger.h>

namespace hydra {

template <typename bit_t, int n_sublat>
GroupActionSublattice<bit_t, n_sublat>::GroupActionSublattice(
    PermutationGroup const &permutation_group)
    : n_sites_(permutation_group.n_sites()),
      n_symmetries_(permutation_group.n_symmetries()),
      permutation_group_(permutation_group),
      n_sites_sublat_(n_sites_ / n_sublat),
      size_tables_(pow(2, n_sites_sublat_)),
      sublat_mask_(((half_bit_t)1 << n_sites_sublat_) - 1) {

  // Check if permutation group is sublattice stable
  if (!symmetries::is_sublattice_stable(n_sublat, permutation_group)) {
    Log.err("Error creating GroupActionSublattice with {} sublattices: "
            "permutation group is not {}-sublattice stable!",
            n_sublat, n_sublat);
  }

  int n_trailing = (n_sublat - 1) * n_sites_sublat_;
  auto group_action = GroupAction(permutation_group_);

  for (int sublat = 0; sublat < n_sublat; ++sublat) {
    sublat_shift_[sublat] = sublat * n_sites_sublat_;

    auto sublat_permutations = symmetries::sublattice_permutations(
        n_sublat, sublat, permutation_group_);
    auto &reps = reps_[sublat];
    auto &rep_syms = rep_syms_[sublat];
    auto &rep_syms_array = rep_syms_array_[sublat];
    reps.resize(size_tables_);
    rep_syms.resize(size_tables_);
    sym_action_[sublat].resize(size_tables_ * n_symmetries_);

    // Compute all representative and representative symmetries
    std::vector<idx_t> rep_syms_begin(size_tables_);
    std::vector<idx_t> rep_syms_end(size_tables_);
    ;
    for (half_bit_t bits = 0; bits < (bit_t)size_tables_; ++bits) {

      idx_t idx = (idx_t)bits;

      bit_t bits_shifted = (bit_t)bits << sublat * n_sites_sublat_;

      // determine the sublattice representatives ...
      bit_t bits_rep = std::numeric_limits<bit_t>::max();
      for (int sym : sublat_permutations) {
        bit_t bits_translated = group_action.apply(sym, bits_shifted);

        // // Security check: bits are translated to least significant bits
        // if (bits == 0) {
        //   assert(bits_translated == 0);
        // } else {
        //   assert(bits_translated >= (bit_t)1 << n_trailing);
        // }

        if (bits_translated < bits_rep) {
          bits_rep = bits_translated;
        }
      }

      // ... and all the symmetries leading to this representative
      rep_syms_begin[idx] = rep_syms_array.size();
      for (int sym : sublat_permutations) {
        bit_t bits_translated = group_action.apply(sym, bits_shifted);
        if (bits_translated == bits_rep) {
          rep_syms_array.push_back(sym);
        }
      }
      rep_syms_end[idx] = rep_syms_array.size();

      // Register new representative and representative symmetries
      half_bit_t rep = (half_bit_t)(bits_rep >> n_trailing);
      reps[idx] = rep;

      // Determine symmetry action on shifted bits
      for (int sym = 0; sym < permutation_group_.size(); ++sym) {
        sym_action(sublat, sym, bits) = group_action.apply(sym, bits_shifted);
      }
    }
    rep_syms_array.shrink_to_fit();

    for (half_bit_t bits = 0; bits < (bit_t)size_tables_; ++bits) {
      idx_t idx = (idx_t)bits;
      rep_syms[idx] =
          gsl::span<int const>(rep_syms_array.data() + rep_syms_begin[idx],
                               rep_syms_end[idx] - rep_syms_begin[idx]);
    }
  }
}

template <typename bit_t, int n_sublat>
bit_t GroupActionSublattice<bit_t, n_sublat>::apply(int sym,
                                                    bit_t state) const {
  bit_t translated = 0;
  for (int sublat = 0; sublat < n_sublat; ++sublat) {
    half_bit_t substate = (state >> sublat_shift_[sublat]) & sublat_mask_;
    translated |= sym_action(sublat, sym, substate);
  }
  return translated;
}

template <typename bit_t, int n_sublat>
bit_t GroupActionSublattice<bit_t, n_sublat>::representative(
    bit_t state) const {

  // Determine sublattice states representatives
  std::array<half_bit_t, n_sublat> sublat_state;
  std::array<half_bit_t, n_sublat> sublat_rep;
  for (int sublat = 0; sublat < n_sublat; ++sublat) {
    sublat_state[sublat] = (state >> sublat_shift_[sublat]) & sublat_mask_;
    sublat_rep[sublat] = reps_[sublat][sublat_state[sublat]];
  }

  // Determine minimal sublattice representative
  half_bit_t min_rep = sublat_rep[0];
  for (int sublat = 1; sublat < n_sublat; ++sublat) {
    if (sublat_rep[sublat] < min_rep) {
      min_rep = sublat_rep[sublat];
    }
  }

  // Apply all sublattice representative symmetries
  bit_t representative = std::numeric_limits<bit_t>::max();
  for (int sublat = 0; sublat < n_sublat; ++sublat) {

    if (sublat_rep[sublat] == min_rep) {

      for (int sym : rep_syms_[sublat][sublat_state[sublat]]) {
        bit_t candidate = 0;

        // Build the translated state
        for (int sl = 0; sl < n_sublat; ++sl) {
          candidate |= sym_action(sl, sym, sublat_state[sl]);
        }

        if (candidate < representative) {
          representative = candidate;
        }
      }
    }
  }
  return representative;
}

template <typename bit_t, int n_sublat>
bool GroupActionSublattice<bit_t, n_sublat>::operator==(
    GroupActionSublattice const &rhs) const {
  return (n_sites_ == rhs.n_sites_) && (n_symmetries_ == rhs.n_symmetries_) &&
         (permutation_group_ == rhs.permutation_group_);
}
template <typename bit_t, int n_sublat>
bool GroupActionSublattice<bit_t, n_sublat>::operator!=(
    GroupActionSublattice const &rhs) const {
  return !operator==(rhs);
}

template class GroupActionSublattice<uint16_t, 1>;
template class GroupActionSublattice<uint32_t, 1>;
template class GroupActionSublattice<uint64_t, 1>;

template class GroupActionSublattice<uint16_t, 2>;
template class GroupActionSublattice<uint32_t, 2>;
template class GroupActionSublattice<uint64_t, 2>;

template class GroupActionSublattice<uint16_t, 3>;
template class GroupActionSublattice<uint32_t, 3>;
template class GroupActionSublattice<uint64_t, 3>;

template class GroupActionSublattice<uint16_t, 4>;
template class GroupActionSublattice<uint32_t, 4>;
template class GroupActionSublattice<uint64_t, 4>;

template class GroupActionSublattice<uint16_t, 5>;
template class GroupActionSublattice<uint32_t, 5>;
template class GroupActionSublattice<uint64_t, 5>;


} // namespace hydra
