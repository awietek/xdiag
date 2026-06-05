// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <set>
#include <string>

#include <xdiag/algebra/rewrite/rules/rules.hpp>

namespace xdiag::algebra {

// MonomialRule: convert the first non-Matrix, non-Id op in a monomial to its
// Matrix form via op_to_matrix_op.
//
// protected_single_op_types — op types that are left as-is when they are the
//   sole op in a size-1 monomial (used by spinhalf_implementation_algebra to
//   keep "Sz", "S+", "S-", "SzSz", "Exchange", "ScalarChirality" in named
//   form). Pass an empty set for matrix_algebra, where every op becomes a
//   Matrix.
//
// Iterated to fixed point together with combine_matrix_rule so all ops in
// multi-op monomials eventually collapse to one Matrix op.
MonomialRule
convert_to_matrix_rule(std::set<std::string> const &protected_single_op_types);

} // namespace xdiag::algebra
