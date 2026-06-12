// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <set>
#include <string>

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// MonomialRule that expands the compound spinless-fermion operators into the
// elementary generators Cdag/C:
//
//   N{i}        -> Cdag{i} C{i}
//   NN{i,j}     -> Cdag{i} C{i} Cdag{j} C{j}
//   Hop{i,j}    -> -Cdag{i} C{j} - Cdag{j} C{i}
//   HopAsym{i,j}-> -Cdag{i} C{j} + Cdag{j} C{i}
//
// A compound operator listed in `protected_types` is left untouched when it is
// the sole operator of a monomial (size 1), so that named operators with a
// dedicated matrix kernel survive and can be dispatched directly. As soon as
// such an operator appears inside a product (size > 1), it is expanded so the
// whole product can be brought into a normal-ordered Cdag/C string.
//
// Operators NOT in `protected_types` are always expanded, even when alone.
// Hence:
//   - fermion_implementation_algebra passes {Hop, HopAsym, N, NN}, keeping the
//     named operators intact for the matrix implementation.
//   - the symmetry-analysis fermion_algebra passes {}, reducing everything to
//     the elementary Cdag/C generators.
//
// The first expandable operator found (left-to-right) is expanded; repeated
// application to a fixed point expands all of them.
MonomialRule
fermion_protected_expansion_rule(std::set<std::string> const &protected_types,
                                 int64_t nsites);

} // namespace xdiag::algebra
