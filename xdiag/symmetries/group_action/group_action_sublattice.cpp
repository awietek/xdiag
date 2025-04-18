// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "group_action_sublattice.hpp"

#include <limits>

#include <xdiag/symmetries/group_action/group_action.hpp>
#include <xdiag/symmetries/group_action/sublattice_stability.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

template <typename bit_t, int n_sublat>
GroupActionSublattice<bit_t, n_sublat>::GroupActionSublattice(
    PermutationGroup const &group)
    : nsites_(group.nsites()), n_symmetries_(group.size()),
      permutation_group_(group), nsites_sublat_(nsites_ / n_sublat),
      size_tables_((int64_t)1 << nsites_sublat_),
      sublat_mask_(((half_bit_t)1 << nsites_sublat_) - 1),
      representative_syms_(n_symmetries_) {

  // Check if permutation group is sublattice stable
  if (!symmetries::is_sublattice_stable(n_sublat, group)) {
    Log.err("Error creating GroupActionSublattice with {} sublattices: "
            "permutation group is not {}-sublattice stable!",
            n_sublat, n_sublat);
  }

  int64_t n_trailing = (n_sublat - 1) * nsites_sublat_;
  auto group_action = GroupAction(permutation_group_);

  for (int sublat = 0; sublat < n_sublat; ++sublat) {
    sublat_shift_[sublat] = sublat * nsites_sublat_;

    auto sublat_permutations = symmetries::sublattice_permutations(
        n_sublat, sublat, permutation_group_);
    auto &reps = reps_[sublat];
    auto &rep_syms = rep_syms_[sublat];
    auto &rep_syms_array = rep_syms_array_[sublat];
    reps.resize(size_tables_);
    rep_syms.resize(size_tables_);
    sym_action_[sublat].resize(size_tables_ * n_symmetries_);

    // Compute all representative and representative symmetries
    std::vector<int64_t> rep_syms_begin(size_tables_);
    std::vector<int64_t> rep_syms_end(size_tables_);

    for (half_bit_t bits = 0; bits < (bit_t)size_tables_; ++bits) {

      int64_t idx = (int64_t)bits;

      bit_t bits_shifted = (bit_t)bits << sublat * nsites_sublat_;

      // determine the sublattice representatives ...
      bit_t bits_rep = std::numeric_limits<bit_t>::max();
      for (int64_t sym : sublat_permutations) {
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
      for (int64_t sym : sublat_permutations) {
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
      for (int64_t sym = 0; sym < permutation_group_.size(); ++sym) {
        sym_action(sublat, sym, bits) = group_action.apply(sym, bits_shifted);
      }
    }
    rep_syms_array.shrink_to_fit();

    for (half_bit_t bits = 0; bits < (bit_t)size_tables_; ++bits) {
      int64_t idx = (int64_t)bits;
      rep_syms[idx] =
          gsl::span<int64_t const>(rep_syms_array.data() + rep_syms_begin[idx],
                                   rep_syms_end[idx] - rep_syms_begin[idx]);
    }
  }
}

template <typename bit_t, int n_sublat>
bit_t GroupActionSublattice<bit_t, n_sublat>::apply(int64_t sym,
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

      for (int64_t sym : rep_syms_[sublat][sublat_state[sublat]]) {
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
std::pair<bit_t, int64_t>
GroupActionSublattice<bit_t, n_sublat>::representative_sym(bit_t state) const {

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
  int64_t representative_sym = 0;
  for (int sublat = 0; sublat < n_sublat; ++sublat) {

    if (sublat_rep[sublat] == min_rep) {

      for (int64_t sym : rep_syms_[sublat][sublat_state[sublat]]) {
        bit_t candidate = 0;

        // Build the translated state
        for (int64_t sl = 0; sl < n_sublat; ++sl) {
          candidate |= sym_action(sl, sym, sublat_state[sl]);
        }

        if (candidate < representative) {
          representative = candidate;
          representative_sym = sym;
        }
      }
    }
  }
  return {representative, representative_sym};
}

template <typename bit_t, int n_sublat>
std::pair<bit_t, gsl::span<int64_t const>>
GroupActionSublattice<bit_t, n_sublat>::representative_syms(bit_t state) const {

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
  gsl::span<int64_t const>::size_type n_syms = 0;
  for (int sublat = 0; sublat < n_sublat; ++sublat) {

    if (sublat_rep[sublat] == min_rep) {

      for (int64_t sym : rep_syms_[sublat][sublat_state[sublat]]) {
        bit_t candidate = 0;

        // Build the translated state
        for (int sl = 0; sl < n_sublat; ++sl) {
          candidate |= sym_action(sl, sym, sublat_state[sl]);
        }

        if (candidate < representative) {
          representative = candidate;
          n_syms = 1;
          representative_syms_[0] = sym;
        } else if (candidate == representative) {
          representative_syms_[n_syms++] = sym;
        }
      }
    }
  }
  return {representative, {representative_syms_.data(), n_syms}};
}

template <typename bit_t, int n_sublat>
bool GroupActionSublattice<bit_t, n_sublat>::operator==(
    GroupActionSublattice const &rhs) const {
  return (nsites_ == rhs.nsites_) && (n_symmetries_ == rhs.n_symmetries_) &&
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

} // namespace xdiag
