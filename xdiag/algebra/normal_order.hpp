// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// Brings an OpSum into normal order with respect to a given Algebra.
//
// A site-free "HubbardU" is treated as an elementary, already normal-ordered
// operator and is left untouched (no expansion into sum_i Nupdn{i}).
//
// Pipeline:
//   1. Expand compound operators into products of elementary ones (OpRules,
//      to fixed point).
//   2. Apply same-site algebra and sort operators into canonical order
//      (MonomialRules, to fixed point). Canonical order: ascending site index;
//      fermionic swaps accumulate a -1 sign per transposition.
//   3. collect: combine equal monomials, drop zero terms.
//   4. collect Matrix ops: merge size-1 Matrix monomials on the same sites
//      into a single Matrix Op (sum of coeff*matrix).
OpSum normal_order(OpSum const &ops, Algebra const &algebra);

} // namespace xdiag::algebra
