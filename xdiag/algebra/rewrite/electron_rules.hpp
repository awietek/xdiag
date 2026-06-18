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

// Spinful-electron rewrite. The expansion is shared by the symmetry and
// implementation algebras; the two differ only in which operators are kept
// named (Algebra::kept_named) and in the normal-order convention.

// Expand one compound electron operator one step toward the elementary
// Cdagup/Cup/Cdagdn/Cdn generators; nullopt if elementary.
std::optional<OpSum> electron_expand(Op const &op, Algebra const &alg);

// Symmetry normal order: site-major sorting + same-site anticommutation, used
// for the permutation-symmetry analysis.
std::optional<OpSum> electron_simplify_symmetry(Monomial const &mono,
                                                Algebra const &alg);

// Implementation normal order: creation-major (all Cdag left, up before dn,
// then C) + same-site anticommutation, expected by the matrix kernels.
std::optional<OpSum> electron_simplify(Monomial const &mono, Algebra const &alg);

} // namespace xdiag::algebra
