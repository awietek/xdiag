// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "spinhalf_sdots_rule.hpp"

#include <vector>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

MonomialRule spinhalf_sdots_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    if (mono.size() != 1 || mono[0].type() != "SdotS") {
      return std::nullopt;
    }
    int64_t i = mono[0][0], j = mono[0][1];
    return OpSum(Op("Exchange", std::vector<int64_t>{i, j})) +
           OpSum(Op("SzSz", std::vector<int64_t>{i, j}));
  };
}

} // namespace xdiag::algebra
