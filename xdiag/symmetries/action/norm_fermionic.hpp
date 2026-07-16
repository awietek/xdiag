// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>

namespace xdiag::symmetries {

// Fermionic counterpart of norm() (symmetries/action/norm.hpp): the norm of a
// state in its symmetry orbit, with every stabilizer character weighted by the
// fermi sign of the permutation acting on the state. Fermi signs can reduce or
// fully cancel the weight, so some bosonically-allowed states acquire zero norm
// and are excluded from the symmetric basis.
template <typename bit_t, typename coeff_t, typename action_t>
double norm_fermionic(bit_t state, action_t const &action,
                      arma::Col<coeff_t> const &characters);

} // namespace xdiag::symmetries
