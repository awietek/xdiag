// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermion_algebra.hpp"

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/fermion_expansion_rules.hpp>
#include <xdiag/algebra/rewrite/rules/fermion_normal_order_rule.hpp>
#include <xdiag/algebra/rewrite/rules/fermion_same_site_rule.hpp>
#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>

namespace xdiag::algebra {

Algebra fermion_algebra(int64_t nsites) {
  std::set<std::string> fermionic{"C", "Cdag"};

  // Empty protected set: every compound operator is always expanded into the
  // elementary Cdag/C generators, even when it is the sole operator of a
  // monomial.
  static const std::set<std::string> protected_types = {};

  std::vector<MonomialRule> algebra_rules_vec{
      id_absorption_rule(),
      fermion_protected_expansion_rule(protected_types),
      fermion_same_site_rule(),
      fermion_normal_order_rule(),
  };

  return Algebra{
      .name = "fermion_algebra",
      .nsites = nsites,
      .d = 2,
      .elementary_types = {"C", "Cdag"},
      .fermionic_types = fermionic,
      .allowed_types = {"C", "Cdag", "Hop", "HopAsym", "Id", "N", "NN"},
      .expansion_rules = {},
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
