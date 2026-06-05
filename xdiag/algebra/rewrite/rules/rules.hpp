// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <functional>
#include <optional>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// A rule that attempts to replace a single Op with an OpSum.
// Returns nullopt if the rule does not apply to the given Op.
using OpRule = std::function<std::optional<OpSum>(Op const &)>;

// A rule that attempts to replace a Monomial (or a sub-sequence of it)
// with an OpSum. Returns nullopt if the rule does not apply.
// The returned OpSum represents the replacement for the matched sub-expression;
// any unmatched prefix/suffix is prepended/appended automatically when
// using pattern_rule. For whole-monomial rules, the returned OpSum replaces
// the entire monomial directly.
using MonomialRule = std::function<std::optional<OpSum>(Monomial const &)>;

} // namespace xdiag::algebra
