// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>

namespace xdiag::symmetries {

// Norm of a state in its symmetry orbit.
// action_t must provide apply(sym, state) and size().
// Explicitly instantiated in norm.cpp for SitePermutation and
// SitePermutationSublattice.
template <typename bit_t, typename coeff_t, typename action_t>
double norm(bit_t state, action_t const &action,
            arma::Col<coeff_t> const &characters);

} // namespace xdiag::symmetries
