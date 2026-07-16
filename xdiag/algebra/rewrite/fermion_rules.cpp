// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "fermion_rules.hpp"

#include <string>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/utils/replace_pair.hpp>
#include <xdiag/algebra/utils/swap_pair.hpp>

namespace xdiag::algebra {

std::optional<OpSum> fermion_expand(Op const &op, Algebra const &alg) {
  std::string t = op.type();

  // Degenerate same-site two-site operators (n in {0,1}): reduce so a kept
  // kernel never sees a same-site bond.
  if (op.hassites() && (op.sites().size() == 2) && (op[0] == op[1])) {
    int64_t s = op[0];
    if (t == "Hop") {
      return -2.0 * Op("N", s);
    }
    if (t == "NN") {
      return OpSum(Op("N", s));
    }
    if (t == "HopAsym") {
      return OpSum{};
    }
  }

  if (t == "N") {
    int64_t i = op[0];
    return OpSum(Op("Cdag", i) * Op("C", i));
  }
  if (t == "TotalN") {
    OpSum r;
    for (int64_t i = 0; i < alg.nsites; ++i) {
      r += Op("N", i);
    }
    return r;
  }
  if (t == "NN") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Cdag", i) * Op("C", i) * Op("Cdag", j) * Op("C", j));
  }
  if (t == "Hop") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdag", i) * Op("C", j));
    r += -1.0 * (Op("Cdag", j) * Op("C", i));
    return r;
  }
  if (t == "HopAsym") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdag", i) * Op("C", j));
    r += +1.0 * (Op("Cdag", j) * Op("C", i));
    return r;
  }
  return std::nullopt;
}

namespace {
bool is_cdagc(std::string const &t) { return (t == "Cdag") || (t == "C"); }
} // namespace

std::optional<OpSum> fermion_simplify(Monomial const &mono, Algebra const &alg) {
  (void)alg;
  int64_t n = mono.size();

  // Pass 1: same-site canonical anticommutation relations.
  for (int64_t k = 0; k + 1 < n; ++k) {
    Op const &a = mono[k];
    Op const &b = mono[k + 1];
    if (a.size() != 1 || b.size() != 1 || a[0] != b[0]) {
      continue;
    }
    std::string ta = a.type(), tb = b.type();
    if (!is_cdagc(ta) || !is_cdagc(tb)) {
      continue;
    }
    int64_t i = a[0];
    if ((ta == "Cdag" && tb == "Cdag") || (ta == "C" && tb == "C")) {
      return replace_pair(mono, k, OpSum{}); // nilpotent
    }
    if (ta == "C" && tb == "Cdag") {
      // {C, Cdag} = Id  =>  C*Cdag = Id - Cdag*C
      return replace_pair(mono, k,
                          OpSum(Op("Id")) - OpSum(Op("Cdag", i) * Op("C", i)));
    }
    // Cdag*C is already canonical.
  }

  // Pass 2: sort different sites into (Cdag before C, ascending site) order.
  for (int64_t k = 0; k + 1 < n; ++k) {
    Op const &a = mono[k];
    Op const &b = mono[k + 1];
    if (a.size() != 1 || b.size() != 1 || a[0] == b[0]) {
      continue;
    }
    std::string ta = a.type(), tb = b.type();
    if (!is_cdagc(ta) || !is_cdagc(tb)) {
      continue;
    }
    int64_t key_a = (ta == "Cdag") ? 0 : 1;
    int64_t key_b = (tb == "Cdag") ? 0 : 1;
    bool out_of_order =
        (key_a > key_b) || ((key_a == key_b) && (a[0] > b[0]));
    if (out_of_order) {
      return swap_pair(mono, k, -1.0); // both fermionic: anticommute
    }
  }
  return std::nullopt;
}

} // namespace xdiag::algebra
