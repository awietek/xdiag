// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "electron_normal_order_rule.hpp"

#include <string>

#include <xdiag/algebra/utils/replace_pair.hpp>
#include <xdiag/algebra/utils/swap_pair.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

namespace {
// Classifies an elementary electron operator type. Returns false for any other
// type. sector: 0 = up, 1 = dn. creation: 0 = Cdag, 1 = C.
bool classify(std::string const &type, int &sector, int &creation) {
  if (type == "Cdagup") {
    sector = 0;
    creation = 0;
  } else if (type == "Cup") {
    sector = 0;
    creation = 1;
  } else if (type == "Cdagdn") {
    sector = 1;
    creation = 0;
  } else if (type == "Cdn") {
    sector = 1;
    creation = 1;
  } else {
    return false;
  }
  return true;
}
} // namespace

MonomialRule electron_normal_order_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    int64_t n = mono.size();
    for (int64_t k = 0; k + 1 < n; ++k) {
      Op const &a = mono[k];
      Op const &b = mono[k + 1];
      if (!a.hassites() || !b.hassites()) {
        continue;
      }
      if (a.size() != 1 || b.size() != 1) {
        continue;
      }
      int sec_a, cre_a, sec_b, cre_b;
      if (!classify(a.type(), sec_a, cre_a) ||
          !classify(b.type(), sec_b, cre_b)) {
        continue; // only elementary electron operators participate
      }
      int64_t sa = a[0], sb = b[0];

      if ((sec_a == sec_b) && (sa == sb)) {
        // Same (sector, site): canonical anticommutation relations.
        if (cre_a == cre_b) {
          // Cdag*Cdag or C*C: nilpotent.
          return replace_pair(mono, k, OpSum{});
        } else if ((cre_a == 0) && (cre_b == 1)) {
          continue; // Cdag*C: number operator, already in order.
        } else {
          // C*Cdag = Id - Cdag*C.
          std::string cdag = (sec_a == 0) ? "Cdagup" : "Cdagdn";
          std::string cann = (sec_a == 0) ? "Cup" : "Cdn";
          OpSum repl = OpSum(Op("Id")) - OpSum(Op(cdag, sa) * Op(cann, sa));
          return replace_pair(mono, k, repl);
        }
      } else {
        // Different modes: anticommute. Canonical key (sector, creation, site).
        bool out_of_order = (sec_a > sec_b) ||
                            ((sec_a == sec_b) && (cre_a > cre_b)) ||
                            ((sec_a == sec_b) && (cre_a == cre_b) && (sa > sb));
        if (out_of_order) {
          return swap_pair(mono, k, -1.0);
        }
      }
    }
    return std::nullopt;
  };
}

} // namespace xdiag::algebra
