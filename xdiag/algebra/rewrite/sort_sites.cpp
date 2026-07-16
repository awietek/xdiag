// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sort_sites.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/utils/swap_pair.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::algebra {

std::optional<OpSum> sort_sites_step(Monomial const &mono, Algebra const &alg) {
  int64_t n = mono.size();
  for (int64_t k = 0; k + 1 < n; ++k) {
    Op const &a = mono[k];
    Op const &b = mono[k + 1];
    if (!a.hassites() || !b.hassites() || a.size() != 1 || b.size() != 1) {
      continue;
    }
    if (a[0] == b[0]) {
      continue; // same site: handled by the block's same-site rule
    }
    if (!(b < a)) {
      continue; // already in order
    }
    bool a_fermi = alg.fermionic_types.count(a.type()) > 0;
    bool b_fermi = alg.fermionic_types.count(b.type()) > 0;
    double sign = (a_fermi && b_fermi) ? -1.0 : 1.0;
    return swap_pair(mono, k, sign);
  }
  return std::nullopt;
}

} // namespace xdiag::algebra
