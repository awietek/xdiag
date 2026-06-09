// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix_algebra.hpp"

#include <set>
#include <string>
#include <vector>

#include <xdiag/algebra/rewrite/rules/boson_expansion_rules.hpp>
#include <xdiag/algebra/rewrite/rules/combine_matrix_rule.hpp>
#include <xdiag/algebra/rewrite/rules/convert_to_matrix_rule.hpp>
#include <xdiag/algebra/rewrite/rules/id_absorption_rule.hpp>
#include <xdiag/algebra/rewrite/rules/sort_matrix_sites_rule.hpp>

namespace xdiag::algebra {

Algebra matrix_algebra(int64_t nsites, int64_t d) {
  std::set<std::string> fermionic{}; // spin operators are bosonic

  // No types are protected: every op (including size-1 Sz, S+, A, Adag, N,
  // etc.) is converted to an explicit Matrix. op_to_matrix_op builds the
  // spin-S / bosonic matrices for the local dimension d, so no OpRules are
  // needed.
  std::vector<MonomialRule> algebra_rules_vec{
      id_absorption_rule(),
      convert_to_matrix_rule({}, d),
      combine_matrix_rule(d),
      sort_matrix_sites_rule(d),
  };

  return Algebra{
      .name = "matrix",
      .nsites = nsites,
      .d = d,
      .elementary_types = {"Matrix"},
      .fermionic_types = fermionic,
      .allowed_types = {"A", "Adag", "Exchange", "ExchangeAsym", "Hop",
                        "HubbardU", "Id", "Matrix", "N", "S+", "S-",
                        "ScalarChirality", "SdotS", "Sx", "Sy", "Sz", "SzSz"},
      .expansion_rules = boson_expansion_rules(nsites),
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::algebra
