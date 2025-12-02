// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>
#include <vector>
#include <xdiag/extern/gsl/span>

namespace xdiag::symmetries {

// determine whether a state is a representative
template <typename bit_t, typename group_action_t>
bool is_representative(bit_t state, group_action_t const &group_action);
  
// Computes the representative (smallest integer value) of "state"
// in orbit given by group_action
template <typename bit_t, typename group_action_t>
bit_t representative(bit_t state, group_action_t const &group_action);

// Computes the representative using only a specified subset of symmetries
template <typename bit_t, class group_action_t>
bit_t representative_subset(bit_t state, group_action_t const &group_action,
                            gsl::span<int64_t const> syms);

// Computes the representative of state AND the symmetry that yields it
template <typename bit_t, class group_action_t>
std::pair<bit_t, int64_t>
representative_sym(bit_t state, group_action_t const &group_action);

// Computes the representative/symmetry using only a specified subset of
// symmetries
template <typename bit_t, class group_action_t>
std::pair<bit_t, int64_t>
representative_sym_subset(bit_t state, group_action_t const &group_action,
                          gsl::span<int64_t const> syms);

// Computes the representative of state and all symmetries that yield it
template <typename bit_t, class group_action_t>
std::pair<bit_t, std::vector<int64_t>>
representative_syms(bit_t state, group_action_t const &group_action);

} // namespace xdiag::symmetries
