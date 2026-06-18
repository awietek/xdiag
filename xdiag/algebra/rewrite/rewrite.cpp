// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "rewrite.hpp"

#include <xdiag/operators/collect.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::algebra {

// One round: apply `step` to every monomial, scaling the replacement by the
// original coefficient. Unchanged monomials are kept. The result is collected.
static OpSum
rewrite_once(OpSum const &ops,
             std::function<std::optional<OpSum>(Monomial const &)> const &step) {
  OpSum result;
  for (auto const &[coeff, mono] : ops) {
    std::optional<OpSum> replacement = step(mono);
    if (replacement) {
      for (auto const &[c_r, m_r] : *replacement) {
        result += OpSum(coeff * c_r, m_r);
      }
    } else {
      result += OpSum(coeff, mono);
    }
  }
  return collect(result);
}

OpSum rewrite_to_fixpoint(
    OpSum const &ops,
    std::function<std::optional<OpSum>(Monomial const &)> const &step,
    int64_t max_iter) try {
  OpSum current = ops;
  for (int64_t iter = 0; iter < max_iter; ++iter) {
    OpSum next = rewrite_once(current, step);
    if (next == current) {
      return next;
    }
    current = std::move(next);
  }
  XDIAG_THROW(fmt::format(
      "rewrite did not reach a fixed point after {} iterations", max_iter));
}
XDIAG_CATCH

} // namespace xdiag::algebra
