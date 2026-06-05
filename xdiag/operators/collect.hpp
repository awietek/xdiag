// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag {

// collect: combine terms with identical monomials by summing their
// coefficients, then drop any terms whose resulting coefficient is zero (within
// tolerance).
//
// As a special case for the "Matrix" type, size-1 "Matrix" monomials acting on
// the same (equally ordered) site vector are fused into a single "Matrix" Op
// whose matrix is the sum of the coefficient-weighted matrices (operator
// linearity); a fused term that sums to the zero matrix is dropped. This is
// well-defined for "Matrix" alone, since a Matrix Op denotes the linear
// operator given by its matrix on those sites, independent of any algebra.
XDIAG_API OpSum collect(OpSum const &ops, double tol = 1e-12);

} // namespace xdiag
