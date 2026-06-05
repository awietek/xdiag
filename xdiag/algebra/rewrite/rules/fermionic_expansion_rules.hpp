// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>
#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// Shared fermionic expansion rules used by electron and tJ algebras: hopping
// (Hop, Hopup, Hopdn), occupation (Nup, Ndn, Ntot, Nupdn), density-density
// pairs, and spin operators in terms of c/c†.
std::vector<OpRule> fermionic_expansion_rules();

} // namespace xdiag::algebra
