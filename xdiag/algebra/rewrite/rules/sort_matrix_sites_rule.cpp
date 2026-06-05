// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sort_matrix_sites_rule.hpp"

#include <algorithm>
#include <vector>

#include <xdiag/algebra/utils/permute_matrix_op.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

MonomialRule sort_matrix_sites_rule(int64_t d) {
  return [d](Monomial const &mono) -> std::optional<OpSum> {
    for (int64_t k = 0; k < mono.size(); ++k) {
      Op const &op = mono[k];
      if (op.type() != "Matrix" || !op.hassites()) {
        continue;
      }
      std::vector<int64_t> const &sites = op.sites();
      if (std::is_sorted(sites.begin(), sites.end())) {
        continue;
      }
      Op permuted = permute_matrix_op(op, d);
      std::vector<Op> ops = mono.ops();
      ops[k] = permuted;
      return OpSum(Monomial(ops));
    }
    return std::nullopt;
  };
}

} // namespace xdiag::algebra
