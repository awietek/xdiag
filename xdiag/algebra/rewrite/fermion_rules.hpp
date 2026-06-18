// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

struct Algebra;

// Spinless-fermion rewrite (used by both the symmetry and the implementation
// algebra -- they share the same Cdag/C normal order and differ only in which
// operators are kept named for a matrix kernel, see Algebra::kept_named).

// Expand one compound operator into simpler ones; nullopt if elementary.
std::optional<OpSum> fermion_expand(Op const &op, Algebra const &alg);

// Same-site canonical anticommutation + sorting (Cdag before C, ascending site)
// of one monomial; nullopt if already normal-ordered.
std::optional<OpSum> fermion_simplify(Monomial const &mono, Algebra const &alg);

} // namespace xdiag::algebra
