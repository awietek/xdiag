// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix_algebra.hpp"

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/combine_matrix_rule.hpp>
#include <xdiag/algebra/rewrite/rules/convert_to_matrix_rule.hpp>
#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>
#include <xdiag/algebra/rewrite/rules/sort_matrix_sites_rule.hpp>

namespace xdiag::algebra {

Algebra matrix_algebra(int64_t d) {
  std::set<std::string> fermionic{}; // spin operators are bosonic

  // No types are protected: every op (including size-1 Sz, S+, etc.) is
  // converted to an explicit Matrix. op_to_matrix_op handles all spin-1/2
  // types directly, so no OpRules are needed.
  std::vector<MonomialRule> algebra_rules_vec{
      id_absorption_rule(),
      convert_to_matrix_rule({}),
      combine_matrix_rule(d),
      sort_matrix_sites_rule(d),
  };

  return Algebra{
      .name = "matrix",
      .d = d,
      .elementary_types = {"Matrix"},
      .fermionic_types = fermionic,
      .allowed_types = {"Exchange", "ExchangeAsym", "Id", "Matrix", "S+", "S-",
                        "ScalarChirality", "SdotS", "Sx", "Sy", "Sz", "SzSz"},
      .expansion_rules = {},
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
