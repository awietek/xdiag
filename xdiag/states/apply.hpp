// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Return the state obtained by applying an operator to the state v
XDIAG_API State apply(Op const &op, State const &v);
XDIAG_API State apply(Monomial const &mono, State const &v);
XDIAG_API State apply(OpSum const &ops, State const &v);

// Apply an operator to v and store the result in the preallocated state w
XDIAG_API void apply(Op const &op, State const &v, State &w);
XDIAG_API void apply(Monomial const &mono, State const &v, State &w);
XDIAG_API void apply(OpSum const &ops, State const &v, State &w);

} // namespace xdiag
