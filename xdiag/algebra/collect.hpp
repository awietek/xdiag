// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag {

// collect: combine terms with identical monomials by summing their coefficients,
// then drop any terms whose resulting coefficient is zero (within tolerance).
// Does NOT reorder monomials — call order() first for canonical ordering.
XDIAG_API OpSum collect(OpSum const &ops, double tol = 1e-12);

} // namespace xdiag
