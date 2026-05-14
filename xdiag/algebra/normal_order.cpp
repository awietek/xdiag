// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "normal_order.hpp"

#include <xdiag/algebra/collect.hpp>
#include <xdiag/algebra/rewrite.hpp>
#include <xdiag/algebra/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::operators {

// Expand site-free HubbardU into sum_{i=0}^{nsites-1} Nupdn{i}
static OpSum expand_hubbard_u(OpSum const &ops, int64_t nsites) {
  OpSum result;
  for (auto const &[coeff, mono] : ops) {
    bool replaced = false;
    for (int64_t k = 0; k < mono.size(); ++k) {
      if (mono[k].type() == "HubbardU" && !mono[k].hassites()) {
        // Build prefix * (sum_i Nupdn{i}) * suffix
        std::vector<Op> pre(mono.ops().begin(), mono.ops().begin() + k);
        std::vector<Op> suf(mono.ops().begin() + k + 1, mono.ops().end());
        Monomial prefix(pre), suffix(suf);

        OpSum nupdn_sum;
        for (int64_t i = 0; i < nsites; ++i) {
          nupdn_sum += Op("Nupdn", i);
        }
        for (auto const &[c2, m2] : nupdn_sum) {
          result += OpSum(coeff * c2, prefix * m2 * suffix);
        }
        replaced = true;
        break; // handle one HubbardU per monomial per pass; rewrite handles the rest
      }
    }
    if (!replaced) {
      result += OpSum(coeff, mono);
    }
  }
  return result;
}

// Validate that all Op types in ops are in the algebra's allowed_types set.
static void check_allowed_types(OpSum const &ops, Algebra const &alg) try {
  if (alg.allowed_types.empty())
    return;
  for (auto const &[coeff, mono] : ops) {
    for (auto const &op : mono) {
      if (alg.allowed_types.count(op.type()) == 0) {
        std::string allowed;
        bool first = true;
        for (auto const &t : alg.allowed_types) {
          if (!first)
            allowed += ", ";
          allowed += t;
          first = false;
        }
        XDIAG_THROW(fmt::format(
            "Op of type \"{}\" is not allowed in the {} algebra.\n"
            "Allowed types: {}",
            op.type(), alg.name, allowed));
      }
    }
  }
}
XDIAG_CATCH

// After collect(), merge size-1 "Matrix" monomials that share the same site
// vector by summing their coefficient-weighted matrices into a single Matrix Op.
static OpSum collect_matrix_ops(OpSum const &ops) {
  using SiteVec = std::vector<int64_t>;
  std::map<SiteVec, Matrix> matrix_sums;
  OpSum non_matrix;

  for (auto const &[coeff, mono] : ops.plain()) {
    if (mono.size() == 1 && mono[0].type() == "Matrix" && mono[0].hassites() &&
        mono[0].hasmatrix()) {
      auto const &sites = mono[0].sites();
      Matrix weighted = mono[0].matrix() * coeff.scalar();
      auto it = matrix_sums.find(sites);
      if (it == matrix_sums.end()) {
        matrix_sums.emplace(sites, std::move(weighted));
      } else {
        it->second += weighted;
      }
    } else {
      non_matrix += OpSum(coeff, mono);
    }
  }

  OpSum result = non_matrix;
  for (auto const &[sites, m] : matrix_sums) {
    result += OpSum(Op("Matrix", sites, m));
  }
  return result;
}

// Validate that all site indices in ops lie in [0, nsites)
static void check_sites_in_range(OpSum const &ops, int64_t nsites) try {
  for (auto const &[coeff, mono] : ops) {
    for (auto const &op : mono) {
      if (op.hassites()) {
        for (int64_t s : op.sites()) {
          if (s < 0 || s >= nsites) {
            XDIAG_THROW(
                fmt::format("Op of type \"{}\" has site index {}, but nsites "
                            "= {}. All sites must be in [0, nsites).",
                            op.type(), s, nsites));
          }
        }
      }
    }
  }
}
XDIAG_CATCH

OpSum normal_order(OpSum const &ops, Algebra const &algebra,
                   int64_t nsites) try {
  // Step 1: validate Op types, expand HubbardU, and validate site indices
  check_allowed_types(ops, algebra);
  OpSum current = expand_hubbard_u(ops.plain(), nsites);
  check_sites_in_range(current, nsites);

  // Step 2: expand compound operators into elementary ones (to fixed point)
  current = rewrite(current, algebra.expansion_rules);

  // Step 3: apply same-site algebra + sort (to fixed point)
  current = rewrite(current, algebra.algebra_rules);

  // Step 4: collect equal monomials
  current = collect(current);

  // Step 5: merge same-sites Matrix ops (sum coeff*matrix per site set)
  return collect_matrix_ops(current);
}
XDIAG_CATCH

} // namespace xdiag::operators
