// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/math/complex.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// computes <v | A | w> expecting real coefficients
XDIAG_API double inner(State const &v, Op const &A, State const &w);
XDIAG_API double inner(State const &v, Monomial const &A, State const &w);
XDIAG_API double inner(State const &v, OpSum const &A, State const &w);

// computes <v | A | w> returning a complex-valed result
XDIAG_API complex innerC(State const &v, Op const &A, State const &w);
XDIAG_API complex innerC(State const &v, Monomial const &A, State const &w);
XDIAG_API complex innerC(State const &v, OpSum const &A, State const &w);

} // namespace xdiag
