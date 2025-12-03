// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <cstdint>
#include <vector>
#include <xdiag/extern/gsl/span>

#pragma once

namespace xdiag::symmetries {

// Computes the norm of a symmetrized state
template <typename bit_t, typename coeff_t, class group_action_t>
double norm(bit_t state, group_action_t const &group_action,
            std::vector<coeff_t> const &characters);

// Computes the norm of a symmetrized state with fermions
template <typename bit_t, typename coeff_t, class group_action_t>
double norm_fermionic(bit_t state, group_action_t const &group_action,
                      std::vector<coeff_t> const &characters);

// Computes the norm of a symmetrized state with up/dn electrons
template <typename bit_t, typename coeff_t, class group_action_t>
double norm_electron(bit_t ups, bit_t dns, group_action_t const &group_action,
                     std::vector<coeff_t> const &characters);

// Computes the norm of a symmetrized state with up/dn electrons (subset of
// syms)
template <typename bit_t, typename coeff_t, class group_action_t>
double norm_electron_subset(bit_t ups, bit_t dns,
                            group_action_t const &group_action,
                            std::vector<coeff_t> const &characters,
                            gsl::span<int64_t const> syms);

} // namespace xdiag::symmetries
