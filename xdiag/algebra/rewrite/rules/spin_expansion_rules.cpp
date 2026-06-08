// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "spin_expansion_rules.hpp"

#include <xdiag/math/complex.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

std::vector<OpRule> spin_expansion_rules() {
  std::vector<OpRule> rules;

  // SdotS{i,j} -> Sz{i}*Sz{j} + 1/2 S+{i}*S-{j} + 1/2 S-{i}*S+{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "SdotS") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += Op("Sz", i) * Op("Sz", j);
    r += 0.5 * (Op("S+", i) * Op("S-", j));
    r += 0.5 * (Op("S-", i) * Op("S+", j));
    return r;
  });

  // SzSz{i,j} -> Sz{i}*Sz{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "SzSz") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Sz", i) * Op("Sz", j));
  });

  // Exchange{i,j} -> 1/2 S+{i}*S-{j} + 1/2 S-{i}*S+{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Exchange") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("S+", i) * Op("S-", j));
    r += 0.5 * (Op("S-", i) * Op("S+", j));
    return r;
  });

  // ExchangeAsym{i,j} -> 1/2 S+{i}*S-{j} - 1/2 S-{i}*S+{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "ExchangeAsym") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("S+", i) * Op("S-", j));
    r -= 0.5 * (Op("S-", i) * Op("S+", j));
    return r;
  });

  // ScalarChirality{i,j,k} ->
  //   (i/2)*[ S+{i}*S-{j}*Sz{k} - S-{i}*S+{j}*Sz{k}
  //         + Sz{i}*S+{j}*S-{k} - Sz{i}*S-{j}*S+{k}
  //         + S-{i}*Sz{j}*S+{k} - S+{i}*Sz{j}*S-{k} ]
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "ScalarChirality") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1], k = op[2];
    complex I(0.0, 1.0);
    OpSum r;
    r += (I * 0.5) * (Op("S+", i) * Op("S-", j) * Op("Sz", k));
    r += (-I * 0.5) * (Op("S-", i) * Op("S+", j) * Op("Sz", k));
    r += (I * 0.5) * (Op("Sz", i) * Op("S+", j) * Op("S-", k));
    r += (-I * 0.5) * (Op("Sz", i) * Op("S-", j) * Op("S+", k));
    r += (I * 0.5) * (Op("S-", i) * Op("Sz", j) * Op("S+", k));
    r += (-I * 0.5) * (Op("S+", i) * Op("Sz", j) * Op("S-", k));
    return r;
  });

  // Sx{i} -> 1/2 S+{i} + 1/2 S-{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sx") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return 0.5 * Op("S+", i) + 0.5 * Op("S-", i);
  });

  // Sy{i} -> -i/2 S+{i} + i/2 S-{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sy") {
      return std::nullopt;
    }
    int64_t i = op[0];
    complex I(0.0, 1.0);
    return (-I * 0.5) * Op("S+", i) + (I * 0.5) * Op("S-", i);
  });

  return rules;
}

} // namespace xdiag::algebra
