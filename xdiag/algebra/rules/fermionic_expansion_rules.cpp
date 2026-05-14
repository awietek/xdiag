// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermionic_expansion_rules.hpp"

#include <xdiag/math/complex.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

std::vector<OpRule> fermionic_expansion_rules() {
  std::vector<OpRule> rules;

  // Hopup{i,j} -> -Cdagup{i}*Cup{j} - Cdagup{j}*Cup{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Hopup") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagup", i) * Op("Cup", j));
    r += -1.0 * (Op("Cdagup", j) * Op("Cup", i));
    return r;
  });

  // Hopdn{i,j} -> -Cdagdn{i}*Cdn{j} - Cdagdn{j}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Hopdn") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagdn", i) * Op("Cdn", j));
    r += -1.0 * (Op("Cdagdn", j) * Op("Cdn", i));
    return r;
  });

  // Hop{i,j} -> -Cdagup{i}*Cup{j} - Cdagup{j}*Cup{i}
  //             -Cdagdn{i}*Cdn{j} - Cdagdn{j}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Hop") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagup", i) * Op("Cup", j));
    r += -1.0 * (Op("Cdagup", j) * Op("Cup", i));
    r += -1.0 * (Op("Cdagdn", i) * Op("Cdn", j));
    r += -1.0 * (Op("Cdagdn", j) * Op("Cdn", i));
    return r;
  });

  // Nup{i} -> Cdagup{i}*Cup{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Nup") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cup", i));
  });

  // Ndn{i} -> Cdagdn{i}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Ndn") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cdn", i));
  });

  // Ntot{i} -> Cdagup{i}*Cup{i} + Cdagdn{i}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Ntot") {
      return std::nullopt;
    }
    int64_t i = op[0];
    OpSum r;
    r += Op("Cdagup", i) * Op("Cup", i);
    r += Op("Cdagdn", i) * Op("Cdn", i);
    return r;
  });

  // Nupdn{i} -> Cdagup{i}*Cup{i}*Cdagdn{i}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Nupdn") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cup", i) * Op("Cdagdn", i) *
                 Op("Cdn", i));
  });

  // NtotNtot{i,j} -> Ntot{i}*Ntot{j}  (further expanded by Ntot rule)
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NtotNtot") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ntot", i) * Op("Ntot", j));
  });

  // NupdnNupdn{i,j} -> Nupdn{i}*Nupdn{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NupdnNupdn") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nupdn", i) * Op("Nupdn", j));
  });

  // NupNdn{i,j} -> Nup{i}*Ndn{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NupNdn") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nup", i) * Op("Ndn", j));
  });

  // NupNup{i,j} -> Nup{i}*Nup{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NupNup") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nup", i) * Op("Nup", j));
  });

  // NdnNdn{i,j} -> Ndn{i}*Ndn{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NdnNdn") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ndn", i) * Op("Ndn", j));
  });

  // NdnNup{i,j} -> Ndn{i}*Nup{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NdnNup") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ndn", i) * Op("Nup", j));
  });

  // Sz{i} -> 1/2 Nup{i} - 1/2 Ndn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sz") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return 0.5 * Op("Nup", i) - 0.5 * Op("Ndn", i);
  });

  // S+{i} -> Cdagup{i} * Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "S+") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cdn", i));
  });

  // S-{i} -> Cdagdn{i} * Cup{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "S-") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cup", i));
  });

  // Sx{i} -> 1/2 S+{i} + 1/2 S-{i}  (further expanded by S+/S- rules above)
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sx") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return 0.5 * Op("S+", i) + 0.5 * Op("S-", i);
  });

  // Sy{i} -> -i/2 S+{i} + i/2 S-{i}  (further expanded by S+/S- rules above)
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

} // namespace xdiag::operators
