// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "spin_algebra.hpp"

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>
#include <xdiag/algebra/rewrite/rules/sort_sites_rule.hpp>
#include <xdiag/algebra/rewrite/rules/spin_expansion_rules.hpp>
#include <xdiag/algebra/rewrite/rules/spinhalf_same_site_rule.hpp>

namespace xdiag::algebra {

Algebra spin_algebra(int64_t nsites) {
  std::set<std::string> fermionic{}; // spin operators commute

  std::vector<MonomialRule> algebra_rules_vec{
      spinhalf_same_site_rule(),
      sort_sites_rule(fermionic),
      id_absorption_rule(),
  };

  return Algebra{
      .name = "spin-1/2",
      .nsites = nsites,
      .d = 2,
      .elementary_types = {"S+", "S-", "Sz"},
      .fermionic_types = fermionic,
      .allowed_types = {"Exchange", "ExchangeAsym", "Id", "S+", "S-",
                        "ScalarChirality", "SdotS", "Sx", "Sy", "Sz", "SzSz"},
      .expansion_rules = spin_expansion_rules(),
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
