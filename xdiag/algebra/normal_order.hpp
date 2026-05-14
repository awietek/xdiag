// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

// Brings an OpSum into normal order with respect to a given Algebra.
//
// nsites: number of lattice sites. Used to:
//   - expand site-free "HubbardU" into sum_{i=0}^{nsites-1} Nupdn{i}
//   - validate that all site indices in ops are in [0, nsites)
//
// Pipeline:
//   1. Expand HubbardU and validate site indices.
//   2. Expand compound operators into products of elementary ones (OpRules,
//      to fixed point).
//   3. Apply same-site algebra and sort operators into canonical order
//      (MonomialRules, to fixed point). Canonical order: ascending site index;
//      fermionic swaps accumulate a -1 sign per transposition.
//   4. collect: combine equal monomials, drop zero terms.
//   5. collect Matrix ops: merge size-1 Matrix monomials on the same sites
//      into a single Matrix Op (sum of coeff*matrix).
OpSum normal_order(OpSum const &ops, Algebra const &algebra, int64_t nsites);

} // namespace xdiag::operators
