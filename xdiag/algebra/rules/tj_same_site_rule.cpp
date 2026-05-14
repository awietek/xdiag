// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj_same_site_rule.hpp"

#include <string>

#include <xdiag/algebra/rules/utils.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

MonomialRule tj_same_site_rule() {
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
      if (a[0] != b[0]) {
        continue;
      }

      int64_t i = a[0];
      std::string ta = a.type(), tb = b.type();

      OpSum repl;
      if (ta == "Cdagup" && tb == "Cdagup") {
        repl = OpSum{};
      } else if (ta == "Cup" && tb == "Cup") {
        repl = OpSum{};
      } else if (ta == "Cdagdn" && tb == "Cdagdn") {
        repl = OpSum{};
      } else if (ta == "Cdn" && tb == "Cdn") {
        repl = OpSum{};
      } else if (ta == "Cup" && tb == "Cdagup") {
        // {Cdagup, Cup} = Id - Ndn  =>  Cup*Cdagup = Id - Cdagup*Cup -
        // Cdagdn*Cdn
        repl = OpSum(Op("Id")) - OpSum(Op("Cdagup", i) * Op("Cup", i)) -
               OpSum(Op("Cdagdn", i) * Op("Cdn", i));
      } else if (ta == "Cdn" && tb == "Cdagdn") {
        // {Cdagdn, Cdn} = Id - Nup  =>  Cdn*Cdagdn = Id - Cdagdn*Cdn -
        // Cdagup*Cup
        repl = OpSum(Op("Id")) - OpSum(Op("Cdagdn", i) * Op("Cdn", i)) -
               OpSum(Op("Cdagup", i) * Op("Cup", i));
      } else if (ta == "Cdn" && tb == "Cdagup") {
        repl = OpSum{}; // 0 (no double occupancy)
      } else if (ta == "Cup" && tb == "Cdagdn") {
        repl = OpSum{}; // 0 (no double occupancy)
      } else if (ta == "Cdagup" && tb == "Cdagdn") {
        repl = -1.0 * (Op("Cdagdn", i) * Op("Cdagup", i)); // Cdagup>Cdagdn
      } else if (ta == "Cup" && tb == "Cdn") {
        repl = -1.0 * (Op("Cdn", i) * Op("Cup", i)); // Cup>Cdn
      } else {
        continue;
      }
      return replace_pair(mono, k, repl);
    }
    return std::nullopt;
  };
}

} // namespace xdiag::operators
