// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <functional>
#include <optional>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// Generic fixed-point driver. `step(mono)` returns the replacement for one
// monomial (without its outer coefficient), or nullopt if that monomial is
// already in final form. rewrite_to_fixpoint applies `step` to every monomial,
// collects equal monomials, and repeats until nothing changes (or throws after
// max_iter rounds). This is the only engine; the per-block expansion and
// simplification steps are plain functions assembled by normal_order().
OpSum rewrite_to_fixpoint(
    OpSum const &ops,
    std::function<std::optional<OpSum>(Monomial const &)> const &step,
    int64_t max_iter = 10000);

} // namespace xdiag::algebra
