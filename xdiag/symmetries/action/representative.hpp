// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>
#include <vector>
#include <xdiag/symmetries/action/site_permutation.hpp>

namespace xdiag::symmetries {

template <typename bit_t>
bit_t representative(bit_t state, SitePermutation const &action);

template <typename bit_t>
std::pair<bit_t, int64_t> representative_sym(bit_t state,
                                             SitePermutation const &action);

template <typename bit_t>
std::pair<bit_t, std::vector<int64_t>>
representative_syms(bit_t state, SitePermutation const &action);

// Representative (smallest value) of `state` using only the given subset of
// symmetries (e.g. an up-stabilizer subgroup for the coupled electron basis).
// `syms` must be non-empty. Works for any bit_t with operator< (initialises
// the running minimum from the first sym rather than numeric_limits).
template <typename bit_t>
bit_t representative_subset(bit_t state, SitePermutation const &action,
                            std::vector<int64_t> const &syms);

// Like representative_subset, but also returns the symmetry that yields it.
template <typename bit_t>
std::pair<bit_t, int64_t>
representative_sym_subset(bit_t state, SitePermutation const &action,
                          std::vector<int64_t> const &syms);

} // namespace xdiag::symmetries
