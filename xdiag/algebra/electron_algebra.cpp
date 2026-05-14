// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "electron_algebra.hpp"

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rules/electron_same_site_rule.hpp>
#include <xdiag/algebra/rules/fermionic_expansion_rules.hpp>
#include <xdiag/algebra/rules/id_absorption_rule.hpp>
#include <xdiag/algebra/rules/sort_sites_rule.hpp>

namespace xdiag::operators {

Algebra electron_algebra() {
  std::set<std::string> fermionic{"Cdagup", "Cup", "Cdagdn", "Cdn"};

  std::vector<MonomialRule> algebra_rules_vec{
      electron_same_site_rule(),
      sort_sites_rule(fermionic),
      id_absorption_rule(),
  };

  return Algebra{
      .name = "electron",
      .elementary_types = {"Cdagup", "Cup", "Cdagdn", "Cdn"},
      .fermionic_types = fermionic,
      .allowed_types = {"Cdagdn", "Cdagup", "Cdn",      "Cup",        "Hop",
                        "Hopdn",  "Hopup",  "HubbardU", "Id",         "Ndn",
                        "NdnNdn", "NdnNup", "Ntot",     "NtotNtot",   "Nup",
                        "NupNdn", "NupNup", "Nupdn",    "NupdnNupdn", "S+",
                        "S-",     "Sx",     "Sy",       "Sz"},
      .expansion_rules = fermionic_expansion_rules(),
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::operators
