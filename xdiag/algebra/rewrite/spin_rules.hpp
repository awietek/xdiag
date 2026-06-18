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

// Spin-1/2 rewrite in terms of S+/S-/Sz (the bosonic spin algebra).

// Expand compound spin operators (SdotS, SzSz, Exchange, ScalarChirality, Sx,
// Sy, TotalSz) into S+/S-/Sz; nullopt if elementary.
std::optional<OpSum> spin_expand(Op const &op, Algebra const &alg);

// Same-site spin-1/2 products + sorting (spin operators commute, no sign);
// nullopt if already normal-ordered.
std::optional<OpSum> spin_simplify(Monomial const &mono, Algebra const &alg);

} // namespace xdiag::algebra
