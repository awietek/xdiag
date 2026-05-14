// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "collect.hpp"

#include <map>
#include <xdiag/math/scalar.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

OpSum collect(OpSum const &ops, double tol) try {
  // Group terms by monomial, summing scalar coefficients.
  // std::map uses Monomial::operator< so the output is in canonical order.
  std::map<Monomial, Scalar> sums;
  for (auto const &[coeff, mono] : ops.plain()) {
    sums[mono] += coeff.scalar();
  }

  OpSum result;
  for (auto const &[mono, c] : sums) {
    if (abs(c) > tol) {
      result += c * mono;
    }
  }
  return result;
}
XDIAG_CATCH

} // namespace xdiag
