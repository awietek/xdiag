// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "combine_matrix_rule.hpp"

#include <xdiag/algebra/utils/combine_matrix_ops.hpp>
#include <xdiag/algebra/utils/replace_pair.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

MonomialRule combine_matrix_rule(int64_t d) {
  return [d](Monomial const &mono) -> std::optional<OpSum> {
    for (int64_t k = 0; k + 1 < mono.size(); ++k) {
      if (mono[k].type() == "Matrix" && mono[k].hassites() &&
          mono[k + 1].type() == "Matrix" && mono[k + 1].hassites()) {
        Op combined = combine_matrix_ops({mono[k], mono[k + 1]}, d);
        return replace_pair(mono, k, OpSum(combined));
      }
    }
    return std::nullopt;
  };
}

} // namespace xdiag::algebra
