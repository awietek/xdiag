// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "electron_protected_expansion_rule.hpp"

#include <optional>
#include <string>
#include <vector>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

namespace {
// Expands one named composite electron operator one step toward the elementary
// generators Cdagup/Cup/Cdagdn/Cdn. Composites built from number operators
// expand to products of Nup/Ndn/Nupdn, which a subsequent pass expands further.
// Returns nullopt for the elementary operators and anything else.
std::optional<OpSum> expand_electron_op(Op const &op, int64_t nsites) {
  std::string type = op.type();

  // Local number operators -> elementary.
  if (type == "Nup") {
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cup", i));
  }
  if (type == "Ndn") {
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cdn", i));
  }
  if (type == "Nupdn") {
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cup", i) * Op("Cdagdn", i) *
                Op("Cdn", i));
  }

  // Hoppings -> elementary.
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

  // Two-site number composites -> products of local number operators.
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
  if (type == "NupdnNupdn") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nupdn", i) * Op("Nupdn", j));
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

  // Site-free interaction.
  if (type == "HubbardU") {
    OpSum r;
    for (int64_t i = 0; i < nsites; ++i) {
      r += Op("Nupdn", i);
    }
    return r;
  }

  return std::nullopt;
}
} // namespace

MonomialRule
electron_protected_expansion_rule(std::set<std::string> const &protected_types,
                                  int64_t nsites) {
  return [protected_types,
          nsites](Monomial const &mono) -> std::optional<OpSum> {
    int64_t n = mono.size();
    for (int64_t k = 0; k < n; ++k) {
      Op const &op = mono[k];

      // Protected named operators survive when they are the sole operator.
      if (n == 1 && protected_types.count(op.type()) > 0) {
        continue;
      }

      std::optional<OpSum> expansion = expand_electron_op(op, nsites);
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
