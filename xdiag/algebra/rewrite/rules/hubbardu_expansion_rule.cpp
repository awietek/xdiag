// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "hubbardu_expansion_rule.hpp"

#include <optional>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

OpRule hubbardu_expansion_rule(int64_t nsites) {
  return [nsites](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "HubbardU") {
      return std::nullopt;
    }
    // HubbardU -> sum_i Nupdn{i}  (Nupdn is the local double-occupancy)
    OpSum r;
    for (int64_t i = 0; i < nsites; ++i) {
      r += Op("Nupdn", i);
    }
    return r;
  };
}

} // namespace xdiag::algebra
