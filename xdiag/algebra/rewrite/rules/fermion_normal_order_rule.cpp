// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermion_normal_order_rule.hpp"

#include <string>

#include <xdiag/algebra/utils/swap_pair.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

MonomialRule fermion_normal_order_rule() {
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
      std::string ta = a.type(), tb = b.type();
      bool a_is_op = (ta == "Cdag") || (ta == "C");
      bool b_is_op = (tb == "Cdag") || (tb == "C");
      if (!a_is_op || !b_is_op) {
        continue; // only Cdag/C participate in fermionic sorting
      }
      if (a[0] == b[0]) {
        continue; // same site: handled by fermion_same_site_rule
      }

      // Canonical order key: Cdag (0) before C (1), then ascending site.
      int64_t key_a = (ta == "Cdag") ? 0 : 1;
      int64_t key_b = (tb == "Cdag") ? 0 : 1;
      bool out_of_order =
          (key_a > key_b) || ((key_a == key_b) && (a[0] > b[0]));
      if (!out_of_order) {
        continue;
      }

      // Both operators are fermionic and on different sites: anticommute.
      return swap_pair(mono, k, -1.0);
    }
    return std::nullopt;
  };
}

} // namespace xdiag::algebra
