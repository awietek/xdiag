// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermion_expansion_rules.hpp"

#include <optional>
#include <string>
#include <vector>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// Expansion of a single compound operator into elementary Cdag/C.
// Returns nullopt for operators that are already elementary (Cdag, C, Id) or
// otherwise not expandable.
static std::optional<OpSum> expand_fermion_op(Op const &op, int64_t nsites) {
  std::string type = op.type();

  // N{i} -> Cdag{i} C{i}
  if (type == "N") {
    int64_t i = op[0];
    return OpSum(Op("Cdag", i) * Op("C", i));
  }

  // TotalN -> sum_i N{i}  (site-free; sums over all `nsites` sites). Each N{i}
  // is expanded to Cdag{i} C{i} by a subsequent pass.
  if (type == "TotalN") {
    OpSum r;
    for (int64_t i = 0; i < nsites; ++i) {
      r += Op("N", i);
    }
    return r;
  }

  // NN{i,j} -> Cdag{i} C{i} Cdag{j} C{j}
  if (type == "NN") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Cdag", i) * Op("C", i) * Op("Cdag", j) * Op("C", j));
  }

  // Hop{i,j} -> -Cdag{i} C{j} - Cdag{j} C{i}
  if (type == "Hop") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdag", i) * Op("C", j));
    r += -1.0 * (Op("Cdag", j) * Op("C", i));
    return r;
  }

  // HopAsym{i,j} -> -Cdag{i} C{j} + Cdag{j} C{i}
  if (type == "HopAsym") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdag", i) * Op("C", j));
    r += +1.0 * (Op("Cdag", j) * Op("C", i));
    return r;
  }

  return std::nullopt;
}

MonomialRule
fermion_protected_expansion_rule(std::set<std::string> const &protected_types,
                                 int64_t nsites) {
  return [protected_types, nsites](Monomial const &mono) -> std::optional<OpSum> {
    int64_t n = mono.size();
    for (int64_t k = 0; k < n; ++k) {
      Op const &op = mono[k];

      // protected named operators survive when they are the sole operator
      if (n == 1 && protected_types.count(op.type()) > 0) {
        continue;
      }

      std::optional<OpSum> expansion = expand_fermion_op(op, nsites);
      if (!expansion) {
        continue;
      }

      // sandwich the expansion between the unmatched prefix and suffix
      std::vector<Op> pre_ops, suf_ops;
      pre_ops.reserve(k);
      suf_ops.reserve(n - k - 1);
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
  };
}

} // namespace xdiag::algebra
