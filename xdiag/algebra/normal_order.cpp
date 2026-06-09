// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "normal_order.hpp"

#include <xdiag/algebra/rewrite/rewrite.hpp>
#include <xdiag/algebra/utils/check_allowed_types.hpp>
#include <xdiag/operators/collect.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::algebra {

OpSum normal_order(OpSum const &ops, Algebra const &algebra) try {
  // Step 1: validate Op types and site ranges
  check_allowed_types(ops, algebra);
  check_sites_in_range(ops, algebra);
  OpSum current = ops.plain();

  // Step 2: expand compound operators into elementary ones (to fixed point)
  current = rewrite(current, algebra.expansion_rules);

  // Step 3: apply same-site algebra + sort (to fixed point)
  current = rewrite(current, algebra.algebra_rules);

  // Step 4: collect equal monomials (also fuses same-site Matrix ops)
  return collect(current);
}
XDIAG_CATCH

} // namespace xdiag::algebra
