// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sort_matrix_sites_rule.hpp"

#include <algorithm>
#include <vector>

#include <xdiag/algebra/permute_matrix_op.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

OpRule sort_matrix_sites_rule(int64_t d) {
  return [d](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Matrix" || !op.hassites()) {
      return std::nullopt;
    }
    std::vector<int64_t> const &sites = op.sites();
    if (std::is_sorted(sites.begin(), sites.end())) {
      return std::nullopt;
    }
    return OpSum(permute_matrix_op(op, d));
  };
}

} // namespace xdiag::operators
