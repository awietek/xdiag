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

// Matrix-path rewrite: operators are turned into explicit local "Matrix"
// operators and fused. Used by the boson / matrix algebra and by the
// spin-1/2 implementation algebra.

// Boson convenience expansions (Hop, HubbardU, TotalN) into A/Adag/N; nullopt
// otherwise (spin operators are left for matrix conversion).
std::optional<OpSum> boson_expand(Op const &op, Algebra const &alg);

// Spin-1/2 implementation expansions (SdotS -> Exchange + SzSz, TotalSz ->
// sum Sz); nullopt otherwise (the kept spin operators have matrix kernels).
std::optional<OpSum> spinhalf_expand(Op const &op, Algebra const &alg);

// Convert each remaining named operator to a local Matrix operator, fuse
// adjacent Matrix operators, and sort their sites. Operators whose type is in
// Algebra::kept_named are left named when they are the sole operator.
std::optional<OpSum> matrix_simplify(Monomial const &mono, Algebra const &alg);

} // namespace xdiag::algebra
