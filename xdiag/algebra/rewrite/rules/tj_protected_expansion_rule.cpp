// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj_protected_expansion_rule.hpp"

#include <optional>
#include <string>
#include <vector>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

namespace {
// Expand one named composite tJ operator one step toward the elementary
// generators Cdagup/Cup/Cdagdn/Cdn. Number composites expand to products of
// Nup/Ndn, which a subsequent pass expands further. Returns nullopt for the
// elementary operators (and anything else).
std::optional<OpSum> tj_expand_op(Op const &op) {
  std::string type = op.type();

  if (type == "Nup") {
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cup", i));
  }
  if (type == "Ndn") {
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cdn", i));
  }
  if (type == "Hopup") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagup", i) * Op("Cup", j));
    r += -1.0 * (Op("Cdagup", j) * Op("Cup", i));
    return r;
  }
  if (type == "Hopdn") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagdn", i) * Op("Cdn", j));
    r += -1.0 * (Op("Cdagdn", j) * Op("Cdn", i));
    return r;
  }
  if (type == "HopupAsym") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagup", i) * Op("Cup", j));
    r += +1.0 * (Op("Cdagup", j) * Op("Cup", i));
    return r;
  }
  if (type == "HopdnAsym") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagdn", i) * Op("Cdn", j));
    r += +1.0 * (Op("Cdagdn", j) * Op("Cdn", i));
    return r;
  }
  if (type == "NupNup") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nup", i) * Op("Nup", j));
  }
  if (type == "NdnNdn") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ndn", i) * Op("Ndn", j));
  }
  if (type == "NupNdn") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nup", i) * Op("Ndn", j));
  }
  if (type == "NdnNup") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ndn", i) * Op("Nup", j));
  }
  if (type == "NtotNtot") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += Op("Nup", i) * Op("Nup", j);
    r += Op("Nup", i) * Op("Ndn", j);
    r += Op("Ndn", i) * Op("Nup", j);
    r += Op("Ndn", i) * Op("Ndn", j);
    return r;
  }
  if (type == "SzSz") {
    // Sz_i Sz_j = 1/4 (Nup_i - Ndn_i)(Nup_j - Ndn_j)
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.25 * (Op("Nup", i) * Op("Nup", j));
    r += -0.25 * (Op("Nup", i) * Op("Ndn", j));
    r += -0.25 * (Op("Ndn", i) * Op("Nup", j));
    r += 0.25 * (Op("Ndn", i) * Op("Ndn", j));
    return r;
  }
  if (type == "tJSzSz") {
    // tJSzSz = SzSz - 1/4 NtotNtot = -1/2 (Nup_i Ndn_j + Ndn_i Nup_j)
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -0.5 * (Op("Nup", i) * Op("Ndn", j));
    r += -0.5 * (Op("Ndn", i) * Op("Nup", j));
    return r;
  }
  if (type == "Exchange") {
    // 1/2 (S+_i S-_j + S-_i S+_j)
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("Cdagup", i) * Op("Cdn", i) * Op("Cdagdn", j) * Op("Cup", j));
    r += 0.5 * (Op("Cdagdn", i) * Op("Cup", i) * Op("Cdagup", j) * Op("Cdn", j));
    return r;
  }
  return std::nullopt;
}
} // namespace

MonomialRule
tj_protected_expansion_rule(std::set<std::string> const &protected_types,
                            int64_t /*nsites*/) {
  return [protected_types](Monomial const &mono) -> std::optional<OpSum> {
    int64_t n = mono.size();
    for (int64_t k = 0; k < n; ++k) {
      Op const &op = mono[k];

      // Protected named operators survive when they are the sole operator.
      if (n == 1 && protected_types.count(op.type()) > 0) {
        continue;
      }

      std::optional<OpSum> expansion = tj_expand_op(op);
      if (!expansion) {
        continue;
      }

      // Sandwich the expansion between the unmatched prefix and suffix.
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
