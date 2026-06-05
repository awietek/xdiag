// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>

#include <xdiag/algebra/rewrite/rules/rules.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// --- Rewriting functions ---

// Apply one round of rewriting to each term in ops:
//   1. Try each MonomialRule on the full monomial — use the first match.
//   2. If no MonomialRule matched, try OpRules left-to-right — use the
//      first match (replaces just that one Op, keeps prefix and suffix).
//   3. If nothing matched, keep the term unchanged.
// The result is collected (equal monomials combined, zero terms dropped).
OpSum rewrite_once(OpSum const &ops, std::vector<MonomialRule> const &mrules,
                   std::vector<OpRule> const &orules);

// Iterate rewrite_once until the OpSum stops changing (fixed point)
// or max_iter is exceeded, in which case XDIAG_THROW is called.
OpSum rewrite(OpSum const &ops, std::vector<MonomialRule> const &mrules,
              std::vector<OpRule> const &orules, int64_t max_iter = 1000);

// Convenience overloads: pass only one kind of rules.
OpSum rewrite_once(OpSum const &ops, std::vector<OpRule> const &orules);
OpSum rewrite_once(OpSum const &ops, std::vector<MonomialRule> const &mrules);
OpSum rewrite(OpSum const &ops, std::vector<OpRule> const &orules,
              int64_t max_iter = 1000);
OpSum rewrite(OpSum const &ops, std::vector<MonomialRule> const &mrules,
              int64_t max_iter = 1000);

} // namespace xdiag::algebra
