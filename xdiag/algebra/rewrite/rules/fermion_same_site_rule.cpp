// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermion_same_site_rule.hpp"

#include <string>

#include <xdiag/algebra/utils/replace_pair.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

MonomialRule fermion_same_site_rule() {
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
      if (ta == "Cdag" && tb == "Cdag") {
        repl = OpSum{};
      } else if (ta == "C" && tb == "C") {
        repl = OpSum{};
      } else if (ta == "C" && tb == "Cdag") {
        // {C, Cdag} = Id  =>  C*Cdag = Id - Cdag*C
        repl = OpSum(Op("Id")) - OpSum(Op("Cdag", i) * Op("C", i));
      } else {
        continue; // Cdag*C is already canonical
      }
      return replace_pair(mono, k, repl);
    }
    return std::nullopt;
  };
}

} // namespace xdiag::algebra
