// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "boson_expansion_rules.hpp"

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

std::vector<OpRule> boson_expansion_rules(int64_t nsites) {
  std::vector<OpRule> rules;

  // Hop{i,j} -> -Adag{i}*A{j} - Adag{j}*A{i}  (minus sign, field convention)
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Hop") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Adag", i) * Op("A", j));
    r += -1.0 * (Op("Adag", j) * Op("A", i));
    return r;
  });

  // HubbardU -> (1/2) sum_i N{i}*(N{i} - 1)
  //           = sum_i [ 1/2 N{i}*N{i} - 1/2 N{i} ]
  // (no sites: sums the on-site interaction over all `nsites` sites)
  rules.push_back([nsites](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "HubbardU") {
      return std::nullopt;
    }
    OpSum r;
    for (int64_t i = 0; i < nsites; ++i) {
      r += 0.5 * (Op("N", i) * Op("N", i));
      r += -0.5 * Op("N", i);
    }
    return r;
  });

  return rules;
}

} // namespace xdiag::algebra
