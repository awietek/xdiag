// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "spin_rules.hpp"

#include <string>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/rewrite/sort_sites.hpp>
#include <xdiag/algebra/utils/replace_pair.hpp>
#include <xdiag/math/complex.hpp>

namespace xdiag::algebra {

std::optional<OpSum> spin_expand(Op const &op, Algebra const &alg) {
  std::string t = op.type();
  if (t == "SdotS") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += Op("Sz", i) * Op("Sz", j);
    r += 0.5 * (Op("S+", i) * Op("S-", j));
    r += 0.5 * (Op("S-", i) * Op("S+", j));
    return r;
  }
  if (t == "SzSz") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Sz", i) * Op("Sz", j));
  }
  if (t == "Exchange") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("S+", i) * Op("S-", j));
    r += 0.5 * (Op("S-", i) * Op("S+", j));
    return r;
  }
  if (t == "ExchangeAsym") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("S+", i) * Op("S-", j));
    r -= 0.5 * (Op("S-", i) * Op("S+", j));
    return r;
  }
  if (t == "ScalarChirality") {
    int64_t i = op[0], j = op[1], k = op[2];
    complex I(0.0, 1.0);
    OpSum r;
    r += (I * 0.5) * (Op("S+", i) * Op("S-", j) * Op("Sz", k));
    r += (-I * 0.5) * (Op("S-", i) * Op("S+", j) * Op("Sz", k));
    r += (I * 0.5) * (Op("Sz", i) * Op("S+", j) * Op("S-", k));
    r += (-I * 0.5) * (Op("Sz", i) * Op("S-", j) * Op("S+", k));
    r += (I * 0.5) * (Op("S-", i) * Op("Sz", j) * Op("S+", k));
    r += (-I * 0.5) * (Op("S+", i) * Op("Sz", j) * Op("S-", k));
    return r;
  }
  if (t == "Sx") {
    int64_t i = op[0];
    return 0.5 * Op("S+", i) + 0.5 * Op("S-", i);
  }
  if (t == "Sy") {
    int64_t i = op[0];
    complex I(0.0, 1.0);
    return (-I * 0.5) * Op("S+", i) + (I * 0.5) * Op("S-", i);
  }
  if (t == "TotalSz") {
    OpSum r;
    for (int64_t i = 0; i < alg.nsites; ++i) {
      r += Op("Sz", i);
    }
    return r;
  }
  return std::nullopt;
}

namespace {
bool is_spin(std::string const &t) {
  return (t == "S+") || (t == "S-") || (t == "Sz");
}
} // namespace

std::optional<OpSum> spin_simplify(Monomial const &mono, Algebra const &alg) {
  (void)alg;
  int64_t n = mono.size();

  // Pass 1: same-site spin-1/2 products.
  for (int64_t k = 0; k + 1 < n; ++k) {
    Op const &a = mono[k];
    Op const &b = mono[k + 1];
    if (a.size() != 1 || b.size() != 1 || a[0] != b[0]) {
      continue;
    }
    std::string ta = a.type(), tb = b.type();
    if (!is_spin(ta) || !is_spin(tb)) {
      continue;
    }
    int64_t i = a[0];
    if (ta == "S+" && tb == "S+") {
      return replace_pair(mono, k, OpSum{});
    }
    if (ta == "S-" && tb == "S-") {
      return replace_pair(mono, k, OpSum{});
    }
    if (ta == "S+" && tb == "Sz") {
      return replace_pair(mono, k, -0.5 * Op("S+", i));
    }
    if (ta == "Sz" && tb == "S+") {
      return replace_pair(mono, k, 0.5 * Op("S+", i));
    }
    if (ta == "S-" && tb == "Sz") {
      return replace_pair(mono, k, 0.5 * Op("S-", i));
    }
    if (ta == "Sz" && tb == "S-") {
      return replace_pair(mono, k, -0.5 * Op("S-", i));
    }
    if (ta == "S+" && tb == "S-") {
      return replace_pair(mono, k, 0.5 * Op("Id") + OpSum(Op("Sz", i)));
    }
    if (ta == "S-" && tb == "S+") {
      return replace_pair(mono, k, 0.5 * Op("Id") - OpSum(Op("Sz", i)));
    }
    if (ta == "Sz" && tb == "Sz") {
      return replace_pair(mono, k, 0.25 * Op("Id"));
    }
  }

  // Pass 2: generic site-major sort (spin operators are bosonic -- no sign).
  return sort_sites_step(mono, alg);
}

} // namespace xdiag::algebra
