// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermion_implementation_algebra.hpp"

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/fermion_expansion_rules.hpp>
#include <xdiag/algebra/rewrite/rules/fermion_normal_order_rule.hpp>
#include <xdiag/algebra/rewrite/rules/fermion_same_site_rule.hpp>
#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>

namespace xdiag::algebra {

Algebra fermion_implementation_algebra(int64_t nsites) {
  std::set<std::string> fermionic{"C", "Cdag"};

  // Named operators that have a dedicated matrix kernel and therefore stay
  // in named form when they appear as a size-1 monomial. As soon as they
  // appear inside a product, they are expanded into elementary Cdag/C.
  static const std::set<std::string> protected_types = {"Hop", "HopAsym", "N",
                                                        "NN"};

  std::vector<MonomialRule> algebra_rules_vec{
      id_absorption_rule(),
      fermion_protected_expansion_rule(protected_types),
      fermion_same_site_rule(),
      fermion_normal_order_rule(),
  };

  return Algebra{
      .name = "fermion_implementation_algebra",
      .nsites = nsites,
      .d = 2,
      .elementary_types = {"C", "Cdag"},
      .fermionic_types = fermionic,
      .allowed_types = {"C", "Cdag", "Hop", "HopAsym", "Id", "N", "NN",
                        "TotalN"},
      .expansion_rules = {},
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
