// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>
#include <xdiag/math/scalar.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {
XDIAG_API bool isapprox(Op const &op1, Op const &op2, double rtol = 1e-12,
                        double atol = 1e-12);
XDIAG_API bool isapprox(OpSum const &ops1, OpSum const &ops2,
                        double rtol = 1e-12, double atol = 1e-12);
std::optional<Scalar> isapprox_multiple(OpSum const &ops1, OpSum const &ops2,
                                        double rtol = 1e-12,
                                        double atol = 1e-12);
} // namespace xdiag
