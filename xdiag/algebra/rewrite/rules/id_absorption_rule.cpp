// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "id_absorption_rule.hpp"

#include <vector>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

MonomialRule id_absorption_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    if (mono.size() <= 1) {
      return std::nullopt;
    }
    for (int64_t k = 0; k < mono.size(); ++k) {
      if (mono[k].type() == "Id") {
        std::vector<Op> ops = mono.ops();
        ops.erase(ops.begin() + k);
        return OpSum(Monomial(ops));
      }
    }
    return std::nullopt;
  };
}

} // namespace xdiag::algebra
