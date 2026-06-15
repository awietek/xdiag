// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj_expansion_rules.hpp"

#include <xdiag/algebra/rewrite/rules/electron_expansion_rules.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

std::vector<OpRule> tj_expansion_rules() {
  std::vector<OpRule> rules = electron_expansion_rules();

  // tJSzSz{i,j} -> (Nup{i}-Ndn{i})/2 * (Nup{j}-Ndn{j})/2
  //  = 1/4*(Nup{i}*Nup{j} - Nup{i}*Ndn{j} - Ndn{i}*Nup{j} + Ndn{i}*Ndn{j})
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "tJSzSz") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.25 * (Op("Nup", i) * Op("Nup", j));
    r += -0.25 * (Op("Nup", i) * Op("Ndn", j));
    r += -0.25 * (Op("Ndn", i) * Op("Nup", j));
    r += 0.25 * (Op("Ndn", i) * Op("Ndn", j));
    return r;
  });

  // tJSdotS{i,j} -> tJSzSz{i,j} + Exchange{i,j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "tJSdotS") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("tJSzSz", {i, j})) + OpSum(Op("Exchange", {i, j}));
  });

  // Exchange{i,j} in tJ: S+{i}*S-{j} + S-{i}*S+{j}
  //   S+{i} = Cdagup{i}*Cdn{i},  S-{i} = Cdagdn{i}*Cup{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Exchange") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    // S+{i}*S-{j} = Cdagup{i}*Cdn{i}*Cdagdn{j}*Cup{j}
    r += Op("Cdagup", i) * Op("Cdn", i) * Op("Cdagdn", j) * Op("Cup", j);
    // S-{i}*S+{j} = Cdagdn{i}*Cup{i}*Cdagup{j}*Cdn{j}
    r += Op("Cdagdn", i) * Op("Cup", i) * Op("Cdagup", j) * Op("Cdn", j);
    return r;
  });

  return rules;
}

} // namespace xdiag::algebra
