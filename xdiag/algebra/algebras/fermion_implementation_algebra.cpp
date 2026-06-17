// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermion_implementation_algebra.hpp"

#include <optional>
#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/fermion_expansion_rules.hpp>
#include <xdiag/algebra/rewrite/rules/fermion_normal_order_rule.hpp>
#include <xdiag/algebra/rewrite/rules/fermion_same_site_rule.hpp>
#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

Algebra fermion_implementation_algebra(int64_t nsites) {
  std::set<std::string> fermionic{"C", "Cdag"};

  // Same-site two-site operators reduce to diagonal single-site operators on the
  // spinless-fermion local space (n in {0,1}): Hop{s,s} = -2 N{s}, NN{s,s} =
  // N{s}, HopAsym{s,s} = 0. Handled in the algebra so the off-diagonal kernel
  // never sees a degenerate same-site bond.
  std::vector<OpRule> expansion_rules;
  expansion_rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (!op.hassites() || op.sites().size() != 2 || op[0] != op[1]) {
      return std::nullopt;
    }
    int64_t s = op[0];
    std::string t = op.type();
    if (t == "Hop") {
      return -2.0 * Op("N", s);
    } else if (t == "NN") {
      return OpSum(Op("N", s));
    } else if (t == "HopAsym") {
      return OpSum();
    }
    return std::nullopt;
  });

  // Named operators that have a dedicated matrix kernel and therefore stay
  // in named form when they appear as a size-1 monomial. As soon as they
  // appear inside a product, they are expanded into elementary Cdag/C (TotalN
  // first expands to sum_i N{i}).
  static const std::set<std::string> protected_types = {"Hop", "HopAsym", "N",
                                                        "NN", "TotalN"};

  std::vector<MonomialRule> algebra_rules_vec{
      id_absorption_rule(),
      fermion_protected_expansion_rule(protected_types, nsites),
      fermion_same_site_rule(),
      fermion_normal_order_rule(),
  };

  return Algebra{
      .name = "fermion_implementation_algebra",
      .nsites = nsites,
      .d = 2,
      .fermionic_types = fermionic,
      .allowed_types = {"C", "Cdag", "Hop", "HopAsym", "Id", "N", "NN",
                        "TotalN"},
      .expansion_rules = expansion_rules,
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
