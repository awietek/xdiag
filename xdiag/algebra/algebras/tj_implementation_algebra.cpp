// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj_implementation_algebra.hpp"

#include <optional>
#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>
#include <xdiag/algebra/rewrite/rules/tj_normal_order_rule.hpp>
#include <xdiag/algebra/rewrite/rules/tj_protected_expansion_rule.hpp>
#include <xdiag/algebra/rewrite/rules/totaln_expansion_rule.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

Algebra tj_implementation_algebra(int64_t nsites) {
  std::set<std::string> fermionic{"Cdagup", "Cup", "Cdagdn", "Cdn"};

  // Named operator types with a dedicated tJ matrix kernel (kept as size-1).
  std::set<std::string> kernel_types{
      "Cdagup", "Cup",    "Cdagdn", "Cdn",      "Hopup",  "Hopdn",
      "Nup",    "Ndn",    "NupNup", "NdnNdn",   "NupNdn", "NdnNup",
      "NtotNtot", "SzSz", "tJSzSz", "Exchange"};

  std::vector<OpRule> expansion_rules;

  // Hop{i,j} -> Hopup{i,j} + Hopdn{i,j}
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Hop") {
      return std::nullopt;
    }
    return OpSum(Op("Hopup", op.sites())) + OpSum(Op("Hopdn", op.sites()));
  });

  // HopAsym{i,j} -> HopupAsym{i,j} + HopdnAsym{i,j}  (antisymmetric / imaginary
  // hopping; the Asym pieces are not kernel types, so the protected-expansion
  // reduces them to Cdag/C strings handled by term_cdagc_string).
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "HopAsym") {
      return std::nullopt;
    }
    return OpSum(Op("HopupAsym", op.sites())) +
           OpSum(Op("HopdnAsym", op.sites()));
  });

  // Same-site two-site operators reduce to diagonal single-site operators on the
  // tJ local space (n in {0,1}, no double occupancy), so the off-diagonal
  // kernels never see a degenerate same-site bond. (Handled in the algebra
  // rather than the kernels; the rewrites are block-specific -- they differ for
  // e.g. the Electron block, which allows double occupancy.)
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (!op.hassites() || op.sites().size() != 2 || op[0] != op[1]) {
      return std::nullopt;
    }
    int64_t s = op[0];
    std::string t = op.type();
    if (t == "Hopup") {
      return -2.0 * Op("Nup", s);
    } else if (t == "Hopdn") {
      return -2.0 * Op("Ndn", s);
    } else if (t == "Exchange") {
      return 0.5 * Op("Ntot", s); // 1/2 (S+S- + S-S+){s,s} = 1/2 n{s}
    } else if (t == "SzSz") {
      return 0.25 * Op("Ntot", s); // Sz{s}^2 = 1/4 n{s}
    } else if (t == "NtotNtot") {
      return OpSum(Op("Ntot", s)); // n{s}^2 = n{s}
    } else if (t == "NupNup") {
      return OpSum(Op("Nup", s));
    } else if (t == "NdnNdn") {
      return OpSum(Op("Ndn", s));
    } else if ((t == "tJSzSz") || (t == "HopupAsym") || (t == "HopdnAsym") ||
               (t == "ExchangeAsym") || (t == "NupNdn") || (t == "NdnNup")) {
      return OpSum(); // identically zero on the tJ local space
    }
    return std::nullopt;
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

  // S+{i} -> Cdagup{i} Cdn{i}
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "S+") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cdn", i));
  });

  // S-{i} -> Cdagdn{i} Cup{i}
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "S-") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cup", i));
  });

  // Sx{i} -> 1/2 S+{i} + 1/2 S-{i}
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sx") {
      return std::nullopt;
    }
    int64_t i = op[0];
    return 0.5 * Op("S+", i) + 0.5 * Op("S-", i);
  });

  // Sy{i} -> -i/2 S+{i} + i/2 S-{i}
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sy") {
      return std::nullopt;
    }
    int64_t i = op[0];
    complex I(0.0, 1.0);
    return (-I * 0.5) * Op("S+", i) + (I * 0.5) * Op("S-", i);
  });

  // SdotS{i,j} -> SzSz{i,j} + Exchange{i,j}   (plain Heisenberg coupling)
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "SdotS") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("SzSz", {i, j})) + OpSum(Op("Exchange", {i, j}));
  });

  // tJSdotS{i,j} -> tJSzSz{i,j} + Exchange{i,j}   (t-J coupling S.S - n n /4)
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "tJSdotS") {
      return std::nullopt;
    }
    int64_t i = op[0], j = op[1];
    return OpSum(Op("tJSzSz", {i, j})) + OpSum(Op("Exchange", {i, j}));
  });

  // ExchangeAsym{i,j} = 1/2 (S+_i S-_j - S-_i S+_j)  (antisymmetric / imaginary
  // exchange; not a kernel type, reduced to Cdag/C strings).
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

  // TotalN -> sum_i Ntot{i}  (Ntot further reduced to Nup + Ndn above)
  expansion_rules.push_back(totaln_expansion_rule(nsites, "Ntot"));

  std::vector<MonomialRule> algebra_rules_vec{
      tj_protected_expansion_rule(kernel_types, nsites),
      tj_normal_order_rule(),
      id_absorption_rule(),
  };

  std::set<std::string> allowed_types = kernel_types;
  for (std::string const &t :
       {"Hop", "HopAsym", "HopupAsym", "HopdnAsym", "Ntot", "Sz", "TotalN",
        "Id", "S+", "S-", "Sx", "Sy", "SdotS", "tJSdotS", "Exchange",
        "ExchangeAsym"}) {
    allowed_types.insert(t);
  }

  return Algebra{
      .name = "tj_implementation",
      .nsites = nsites,
      .d = 3,
      .fermionic_types = fermionic,
      .allowed_types = allowed_types,
      .expansion_rules = expansion_rules,
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
