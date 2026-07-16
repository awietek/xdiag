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

// tJ rewrite (electron space projected to no double occupancy). The expansion
// is shared by the symmetry and implementation algebras; the normal-order
// conventions differ.

// Expand one compound tJ operator one step toward Cdagup/Cup/Cdagdn/Cdn;
// nullopt if elementary.
std::optional<OpSum> tj_expand(Op const &op, Algebra const &alg);

// Symmetry normal order: site-major sorting + projected same-site relations.
std::optional<OpSum> tj_simplify_symmetry(Monomial const &mono,
                                          Algebra const &alg);

// Implementation normal order: creation-major + projected same-site relations
// (uses the tJ-specific (anti)commutators, not the canonical ones).
std::optional<OpSum> tj_simplify(Monomial const &mono, Algebra const &alg);

} // namespace xdiag::algebra
