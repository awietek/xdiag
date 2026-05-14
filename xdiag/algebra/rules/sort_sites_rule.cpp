// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sort_sites_rule.hpp"

#include <xdiag/algebra/rules/utils.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

MonomialRule sort_sites_rule(std::set<std::string> const &fermionic_types) {
  return [fermionic_types](Monomial const &mono) -> std::optional<OpSum> {
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
      if (a[0] == b[0]) {
        continue; // same-site: handled by algebra rules, not here
      }
      if (!(b < a)) {
        continue; // already in order
      }

      bool a_fermi = fermionic_types.count(a.type()) > 0;
      bool b_fermi = fermionic_types.count(b.type()) > 0;
      double sign = (a_fermi && b_fermi) ? -1.0 : 1.0;
      return swap_pair(mono, k, sign);
    }
    return std::nullopt;
  };
}

} // namespace xdiag::operators
