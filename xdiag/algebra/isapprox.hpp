// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

bool isapprox(Op const &op1, Op const &op2, double rtol = 1e-12,
              double atol = 1e-12);

// Brings both OpSums into normal order with respect to the given algebra and
// then compares the normal-ordered forms.
bool isapprox(OpSum const &ops1, OpSum const &ops2,
              algebra::Algebra const &algebra, double rtol = 1e-12,
              double atol = 1e-12);

// Brings both OpSums into normal order with respect to the given algebra and
// returns lambda such that ops1 = lambda * ops2, or std::nullopt if no such
// scalar exists.
std::optional<Scalar> isapprox_multiple(OpSum const &ops1, OpSum const &ops2,
                                        algebra::Algebra const &algebra,
                                        double rtol = 1e-12,
                                        double atol = 1e-12);

} // namespace xdiag
