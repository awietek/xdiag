// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "normal_order.hpp"

#include <optional>
#include <vector>

#include <xdiag/algebra/rewrite/rewrite.hpp>
#include <xdiag/algebra/utils/check_allowed_types.hpp>
#include <xdiag/operators/collect.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::algebra {

// Expand one monomial by one step: find the first operator that should expand
// and substitute its expansion (sandwiched between the surrounding operators).
// A size-1 operator whose type is in kept_named survives -- it has a dedicated
// matrix kernel -- unless it is a degenerate same-site two-site operator (e.g.
// Hop{s,s}), which always reduces.
static std::optional<OpSum> expand_monomial(Monomial const &mono,
                                            Algebra const &alg) {
  if (!alg.expand) {
    return std::nullopt;
  }
  int64_t n = mono.size();
  for (int64_t k = 0; k < n; ++k) {
    Op const &op = mono[k];
    bool degenerate = op.hassites() && (op.sites().size() == 2) &&
                      (op.sites()[0] == op.sites()[1]);
    if ((n == 1) && (alg.kept_named.count(op.type()) > 0) && !degenerate) {
      continue; // keep named for its kernel
    }
    std::optional<OpSum> expansion = alg.expand(op, alg);
    if (!expansion) {
      continue;
    }
    std::vector<Op> pre_ops, suf_ops;
    for (int64_t j = 0; j < k; ++j) {
      pre_ops.push_back(mono[j]);
    }
    for (int64_t j = k + 1; j < n; ++j) {
      suf_ops.push_back(mono[j]);
    }
    Monomial prefix(pre_ops), suffix(suf_ops);
    OpSum result;
    for (auto const &[c_r, m_r] : *expansion) {
      result += OpSum(c_r, prefix * m_r * suffix);
    }
    return result;
  }
  return std::nullopt;
}

// Remove one Id factor from a product (an Id in a size-1 monomial is the
// identity operator itself and is kept). Generic, block-independent.
static std::optional<OpSum> remove_id(Monomial const &mono) {
  if (mono.size() <= 1) {
    return std::nullopt;
  }
  for (int64_t k = 0; k < mono.size(); ++k) {
    if (mono[k].type() == "Id") {
      std::vector<Op> ops = mono.ops();
      ops.erase(ops.begin() + k);
      return OpSum(Monomial(ops));
    }
  }
  return std::nullopt;
}

// Simplify one monomial by one step: first absorb a stray Id, then defer to the
// block's same-site / sorting rule.
static std::optional<OpSum> simplify_monomial(Monomial const &mono,
                                              Algebra const &alg) {
  std::optional<OpSum> id_removed = remove_id(mono);
  if (id_removed) {
    return id_removed;
  }
  if (alg.simplify) {
    return alg.simplify(mono, alg);
  }
  return std::nullopt;
}

OpSum normal_order(OpSum const &ops, Algebra const &algebra) try {
  // Step 1: validate Op types and site ranges.
  check_allowed_types(ops, algebra);
  check_sites_in_range(ops, algebra);
  OpSum current = ops.plain();

  // Step 2: expand compound operators to a fixed point.
  current = rewrite_to_fixpoint(current, [&](Monomial const &mono) {
    return expand_monomial(mono, algebra);
  });

  // Step 3: same-site simplification + sorting to a fixed point.
  current = rewrite_to_fixpoint(current, [&](Monomial const &mono) {
    return simplify_monomial(mono, algebra);
  });

  // Step 4: collect equal monomials (also fuses same-site Matrix ops).
  return collect(current);
}
XDIAG_CATCH

} // namespace xdiag::algebra
