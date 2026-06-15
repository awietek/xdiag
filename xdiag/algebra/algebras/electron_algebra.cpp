// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "electron_algebra.hpp"

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/electron_expansion_rules.hpp>
#include <xdiag/algebra/rewrite/rules/electron_same_site_rule.hpp>
#include <xdiag/algebra/rewrite/rules/hubbardu_expansion_rule.hpp>
#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>
#include <xdiag/algebra/rewrite/rules/sort_sites_rule.hpp>
#include <xdiag/algebra/rewrite/rules/totaln_expansion_rule.hpp>

namespace xdiag::algebra {

Algebra electron_algebra(int64_t nsites) {
  std::set<std::string> fermionic{"Cdagup", "Cup", "Cdagdn", "Cdn"};

  std::vector<MonomialRule> algebra_rules_vec{
      electron_same_site_rule(),
      sort_sites_rule(fermionic),
      id_absorption_rule(),
  };

  std::vector<OpRule> expansion_rules_vec = electron_expansion_rules();
  expansion_rules_vec.push_back(totaln_expansion_rule(nsites, "Ntot"));
  expansion_rules_vec.push_back(hubbardu_expansion_rule(nsites));

  return Algebra{
      .name = "electron",
      .nsites = nsites,
      .d = 4,
      .fermionic_types = fermionic,
      .allowed_types = {"Cdagdn",     "Cdagup",   "Cdn",       "Cup",
                        "Exchange",   "ExchangeAsym", "Hop",   "HopAsym",
                        "Hopdn",      "HopdnAsym", "Hopup",    "HopupAsym",
                        "HubbardU",   "Id",       "Ndn",       "NdnNdn",
                        "NdnNup",     "Ntot",     "NtotNtot",  "Nup",
                        "NupNdn",     "NupNup",   "Nupdn",     "NupdnNupdn",
                        "S+",         "S-",       "SdotS",     "Sx",
                        "Sy",         "Sz",       "SzSz",      "TotalN"},
      .expansion_rules = expansion_rules_vec,
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
