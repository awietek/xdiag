// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj_normal_order_rule.hpp"

#include <string>

#include <xdiag/algebra/utils/replace_pair.hpp>
#include <xdiag/algebra/utils/swap_pair.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

namespace {
// sector: 0 = up, 1 = dn. creation: 0 = Cdag, 1 = C. Returns false otherwise.
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

MonomialRule tj_normal_order_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    int64_t n = mono.size();
    for (int64_t k = 0; k + 1 < n; ++k) {
      Op const &a = mono[k];
      Op const &b = mono[k + 1];
      if (!a.hassites() || !b.hassites() || a.size() != 1 || b.size() != 1) {
        continue;
      }
      int sec_a, cre_a, sec_b, cre_b;
      if (!classify(a.type(), sec_a, cre_a) ||
          !classify(b.type(), sec_b, cre_b)) {
        continue;
      }
      int64_t sa = a[0], sb = b[0];

      if (sa != sb) {
        // Different sites: ordinary fermion anticommutation. Creation-major
        // canonical order -- all creation operators left of all annihilation
        // operators, then up before dn, then ascending site: key (creation,
        // sector, site). (The same-site branch below already produces this
        // order: it keeps Cdag*C pairs S+/S- and forbids the C*Cdag ones.)
        bool out_of_order = (cre_a > cre_b) ||
                            ((cre_a == cre_b) && (sec_a > sec_b)) ||
                            ((cre_a == cre_b) && (sec_a == sec_b) && (sa > sb));
        if (out_of_order) {
          return swap_pair(mono, k, -1.0);
        }
        continue;
      }

      // Same site i: tJ (projected) relations.
      int64_t i = sa;
      if (sec_a == sec_b) {
        if (cre_a == cre_b) {
          return replace_pair(mono, k, OpSum{}); // Cdag*Cdag / C*C nilpotent
        } else if ((cre_a == 0) && (cre_b == 1)) {
          continue; // Cdag*C: number operator, already in order
        } else {
          // C*Cdag: {C,Cdag} = Id - n_other  =>
          // C*Cdag = Id - Cdag*C - Cdag_other*C_other
          std::string cdag = (sec_a == 0) ? "Cdagup" : "Cdagdn";
          std::string cann = (sec_a == 0) ? "Cup" : "Cdn";
          std::string cdag_o = (sec_a == 0) ? "Cdagdn" : "Cdagup";
          std::string cann_o = (sec_a == 0) ? "Cdn" : "Cup";
          OpSum repl = OpSum(Op("Id")) - OpSum(Op(cdag, i) * Op(cann, i)) -
                       OpSum(Op(cdag_o, i) * Op(cann_o, i));
          return replace_pair(mono, k, repl);
        }
      } else {
        // Different sector, same site.
        if (cre_a == cre_b) {
          // Both creation or both annihilation: order up(0) before dn(1).
          if (sec_a > sec_b) {
            return swap_pair(mono, k, -1.0);
          }
          continue;
        } else if ((cre_a == 0) && (cre_b == 1)) {
          // Cdag then C, opposite spin: S+ (Cdagup Cdn) / S- (Cdagdn Cup).
          // Both are kept as-is (canonical tJ form).
          continue;
        } else {
          // C then Cdag, opposite spin: Cup*Cdagdn / Cdn*Cdagup = 0 (a creation
          // onto the spin already blocked by the other spin's projection).
          return replace_pair(mono, k, OpSum{});
        }
      }
    }
    return std::nullopt;
  };
}

} // namespace xdiag::algebra
