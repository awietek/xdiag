// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>
#include <vector>

#include <xdiag/symmetries/operations/fermi_sign.hpp>

namespace xdiag::symmetries {

// Fermionic counterpart of norm() (symmetries/action/norm.hpp): the norm of a
// state in its symmetry orbit, with every stabilizer character weighted by the
// fermi sign of the permutation acting on the state. Fermi signs can reduce or
// fully cancel the weight, so some bosonically-allowed states acquire zero norm
// and are excluded from the symmetric basis.
//
// Header-only template: it is instantiated for every bit_t reached through the
// representative table's fermionic branch, which would otherwise need a long
// explicit-instantiation list.
template <typename bit_t, typename coeff_t, typename action_t>
double norm_fermionic(bit_t state, action_t const &action,
                      arma::Col<coeff_t> const &characters) {
  auto const &group = action.group();
  coeff_t amplitude = 0.0;
  for (int64_t sym = 0; sym < action.size(); ++sym) {
    if (action.apply(sym, state) == state) {
      bool fermi = fermi_bool_of_permutation(state, group[sym]);
      amplitude += fermi ? -characters(sym) : characters(sym);
    }
  }
  return std::sqrt(std::abs(amplitude));
}

} // namespace xdiag::symmetries
