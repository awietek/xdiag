// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "spinhalf_implementation_algebra.hpp"

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/combine_matrix_rule.hpp>
#include <xdiag/algebra/rewrite/rules/convert_to_matrix_rule.hpp>
#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>
#include <xdiag/algebra/rewrite/rules/sort_matrix_sites_rule.hpp>
#include <xdiag/algebra/rewrite/rules/spinhalf_sdots_rule.hpp>

namespace xdiag::algebra {

Algebra spinhalf_implementation_algebra(int64_t nsites) {
  std::set<std::string> fermionic{};

  // Types that remain in named form when they appear as a size-1 monomial.
  // All other allowed types (Sx, Sy, SdotS after sdots_rule, multi-op
  // combinations) are converted to Matrix.
  static const std::set<std::string> protected_types = {
      "Exchange", "ExchangeAsym", "S+", "S-", "ScalarChirality", "Sz", "SzSz"};

  std::vector<MonomialRule> algebra_rules_vec{
      id_absorption_rule(),
      spinhalf_sdots_rule(),
      convert_to_matrix_rule(protected_types, 2),
      combine_matrix_rule(2),
      sort_matrix_sites_rule(2),
  };

  return Algebra{
      .name = "spinhalf_implementation_algebra",
      .nsites = nsites,
      .d = 2,
      .elementary_types = {"Exchange", "ExchangeAsym", "Matrix", "S+", "S-",
                           "ScalarChirality", "Sz", "SzSz"},
      .fermionic_types = fermionic,
      .allowed_types = {"Exchange", "ExchangeAsym", "Id", "Matrix", "S+", "S-",
                        "ScalarChirality", "SdotS", "Sx", "Sy", "Sz", "SzSz"},
      .expansion_rules = {},
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
