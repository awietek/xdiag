// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "electron_implementation_algebra.hpp"

#include <optional>
#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/electron_normal_order_rule.hpp>
#include <xdiag/algebra/rewrite/rules/electron_protected_expansion_rule.hpp>
#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>
#include <xdiag/algebra/rewrite/rules/totaln_expansion_rule.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

Algebra electron_implementation_algebra(int64_t nsites) {
  std::set<std::string> fermionic{"Cdagup", "Cup", "Cdagdn", "Cdn"};

  // Named operator types that have a dedicated matrix kernel
  // (matrices/blocks/electron/terms). They are kept as size-1 named ops.
  std::set<std::string> kernel_types{
      "Cdagup", "Cup",        "Cdagdn",   "Cdn",    "Hopup",
      "Hopdn",  "HopupAsym",  "HopdnAsym", "Nup",   "Ndn",
      "Nupdn",  "NupdnNupdn", "NtotNtot",  "NupNdn", "NupNup",
      "NdnNdn", "NdnNup",     "SzSz",      "HubbardU"};

  // Convenience operators are reduced to the kernel types.
  std::vector<OpRule> expansion_rules;

  // Hop{i,j} -> Hopup{i,j} + Hopdn{i,j}
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Hop") {
      return std::nullopt;
    }
    return OpSum(Op("Hopup", op.sites())) + OpSum(Op("Hopdn", op.sites()));
  });

  // HopAsym{i,j} -> HopupAsym{i,j} + HopdnAsym{i,j}
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "HopAsym") {
      return std::nullopt;
    }
    return OpSum(Op("HopupAsym", op.sites())) +
           OpSum(Op("HopdnAsym", op.sites()));
  });

  // Ntot{i} -> Nup{i} + Ndn{i}
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Ntot") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Nup", i)) + OpSum(Op("Ndn", i));
  });

  // Sz{i} -> 1/2 Nup{i} - 1/2 Ndn{i}
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sz") {
      return std::nullopt;
    }
    int64_t i = op[0];
    OpSum r;
    r += 0.5 * Op("Nup", i);
    r += -0.5 * Op("Ndn", i);
    return r;
  });

  // S+{i} -> Cdagup{i} Cdn{i}   (already "all ups then all dns" order)
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "S+") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cdn", i));
  });

  // S-{i} -> Cdagdn{i} Cup{i}  (the normal-order rule reorders to ups-left)
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "S-") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cup", i));
  });

  // Sx{i} -> 1/2 S+{i} + 1/2 S-{i}  (further reduced by the S+/S- rules)
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sx") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return 0.5 * Op("S+", i) + 0.5 * Op("S-", i);
  });

  // Sy{i} -> -i/2 S+{i} + i/2 S-{i}  (further reduced by the S+/S- rules)
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sy") {
      return std::nullopt;
    }
    int64_t i = op[0];
    complex I(0.0, 1.0);
    return (-I * 0.5) * Op("S+", i) + (I * 0.5) * Op("S-", i);
  });

  // Exchange{i,j} = 1/2 (S+_i S-_j + S-_i S+_j)
  //   = 1/2 Cdagup_i Cdn_i Cdagdn_j Cup_j + 1/2 Cdagdn_i Cup_i Cdagup_j Cdn_j
  // (the normal-order rule then brings each product to ups-left form).
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Exchange") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("Cdagup", i) * Op("Cdn", i) * Op("Cdagdn", j) * Op("Cup", j));
    r += 0.5 * (Op("Cdagdn", i) * Op("Cup", i) * Op("Cdagup", j) * Op("Cdn", j));
    return r;
  });

  // ExchangeAsym{i,j} = 1/2 (S+_i S-_j - S-_i S+_j)
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "ExchangeAsym") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("Cdagup", i) * Op("Cdn", i) * Op("Cdagdn", j) * Op("Cup", j));
    r -= 0.5 * (Op("Cdagdn", i) * Op("Cup", i) * Op("Cdagup", j) * Op("Cdn", j));
    return r;
  });

  // SdotS{i,j} -> SzSz{i,j} + Exchange{i,j}  (SzSz kept as a kernel, Exchange
  // further reduced to a Cdag/C string above)
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "SdotS") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("SzSz", {i, j})) + OpSum(Op("Exchange", {i, j}));
  });

  // TotalN -> sum_i Ntot{i}  (Ntot is further reduced to Nup + Ndn above)
  expansion_rules.push_back(totaln_expansion_rule(nsites, "Ntot"));

  std::vector<MonomialRule> algebra_rules_vec{
      // Composite named operators expand to elementary Cdag/C strings when they
      // appear inside a product (kept named when alone -> dedicated kernel).
      electron_protected_expansion_rule(kernel_types, nsites),
      electron_normal_order_rule(),
      id_absorption_rule(),
  };

  std::set<std::string> allowed_types = kernel_types;
  for (std::string const &t :
       {"Hop", "HopAsym", "Ntot", "Sz", "TotalN", "Id", "S+", "S-", "Sx", "Sy",
        "Exchange", "ExchangeAsym", "SdotS"}) {
    allowed_types.insert(t);
  }

  return Algebra{
      .name = "electron_implementation",
      .nsites = nsites,
      .d = 4,
      .fermionic_types = fermionic,
      .allowed_types = allowed_types,
      .expansion_rules = expansion_rules,
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
