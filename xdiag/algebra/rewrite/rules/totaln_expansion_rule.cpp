// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "totaln_expansion_rule.hpp"

#include <optional>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

OpRule totaln_expansion_rule(int64_t nsites) {
  return [nsites](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "TotalN") {
      return std::nullopt;
    }
    // TotalN -> sum_i N{i}  (no sites: sums over all `nsites` sites)
    OpSum r;
    for (int64_t i = 0; i < nsites; ++i) {
      r += Op("N", i);
    }
    return r;
  };
}

} // namespace xdiag::algebra
