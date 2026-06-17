// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj_algebra.hpp"

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>
#include <xdiag/algebra/rewrite/rules/sort_sites_rule.hpp>
#include <xdiag/algebra/rewrite/rules/tj_expansion_rules.hpp>
#include <xdiag/algebra/rewrite/rules/tj_same_site_rule.hpp>
#include <xdiag/algebra/rewrite/rules/totaln_expansion_rule.hpp>

namespace xdiag::algebra {

Algebra tj_algebra(int64_t nsites) {
  std::set<std::string> fermionic{"Cdagup", "Cup", "Cdagdn", "Cdn"};

  std::vector<MonomialRule> algebra_rules_vec{
      tj_same_site_rule(),
      sort_sites_rule(fermionic),
      id_absorption_rule(),
  };

  // Nupdn / NupdnNupdn are identically zero on the tJ space (no double
  // occupancy), so they are not allowed types. TotalN = sum_i Ntot{i}.
  std::vector<OpRule> expansion_rules = tj_expansion_rules();
  expansion_rules.push_back(totaln_expansion_rule(nsites, "Ntot"));

  return Algebra{
      .name = "tJ",
      .nsites = nsites,
      .d = 3,
      .fermionic_types = fermionic,
      .allowed_types = {"Cdagdn",   "Cdagup",       "Cdn",      "Cup",
                        "Exchange", "ExchangeAsym", "Hop",      "HopAsym",
                        "Hopdn",    "HopdnAsym",    "Hopup",    "HopupAsym",
                        "Id",       "Ndn",          "NdnNdn",   "NdnNup",
                        "Ntot",     "NtotNtot",     "Nup",      "NupNdn",
                        "NupNup",   "S+",           "S-",       "SdotS",
                        "SzSz",     "Sx",           "Sy",       "Sz",
                        "TotalN",   "tJSdotS",      "tJSzSz"},
      .expansion_rules = expansion_rules,
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
