// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "tj_rules.hpp"

#include <string>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/rewrite/sort_sites.hpp>
#include <xdiag/algebra/utils/replace_pair.hpp>
#include <xdiag/algebra/utils/swap_pair.hpp>
#include <xdiag/math/complex.hpp>

namespace xdiag::algebra {

std::optional<OpSum> tj_expand(Op const &op, Algebra const &alg) {
  std::string t = op.type();

  // Degenerate same-site two-site operators (n in {0,1}, no double occupancy).
  if (op.hassites() && (op.sites().size() == 2) && (op[0] == op[1])) {
    int64_t s = op[0];
    if (t == "Hopup") {
      return -2.0 * Op("Nup", s);
    }
    if (t == "Hopdn") {
      return -2.0 * Op("Ndn", s);
    }
    if (t == "Exchange") {
      return 0.5 * Op("Ntot", s); // 1/2 (S+S- + S-S+){s,s} = 1/2 n{s}
    }
    if (t == "SzSz") {
      return 0.25 * Op("Ntot", s); // Sz{s}^2 = 1/4 n{s}
    }
    if (t == "NtotNtot") {
      return OpSum(Op("Ntot", s));
    }
    if (t == "NupNup") {
      return OpSum(Op("Nup", s));
    }
    if (t == "NdnNdn") {
      return OpSum(Op("Ndn", s));
    }
    if ((t == "tJSzSz") || (t == "HopupAsym") || (t == "HopdnAsym") ||
        (t == "ExchangeAsym") || (t == "NupNdn") || (t == "NdnNup")) {
      return OpSum{}; // identically zero on the tJ local space
    }
  }

  if (t == "Hop") {
    return OpSum(Op("Hopup", op.sites())) + OpSum(Op("Hopdn", op.sites()));
  }
  if (t == "HopAsym") {
    return OpSum(Op("HopupAsym", op.sites())) +
           OpSum(Op("HopdnAsym", op.sites()));
  }
  if (t == "Hopup") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagup", i) * Op("Cup", j));
    r += -1.0 * (Op("Cdagup", j) * Op("Cup", i));
    return r;
  }
  if (t == "Hopdn") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagdn", i) * Op("Cdn", j));
    r += -1.0 * (Op("Cdagdn", j) * Op("Cdn", i));
    return r;
  }
  if (t == "HopupAsym") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagup", i) * Op("Cup", j));
    r += +1.0 * (Op("Cdagup", j) * Op("Cup", i));
    return r;
  }
  if (t == "HopdnAsym") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagdn", i) * Op("Cdn", j));
    r += +1.0 * (Op("Cdagdn", j) * Op("Cdn", i));
    return r;
  }

  if (t == "Nup") {
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cup", i));
  }
  if (t == "Ndn") {
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cdn", i));
  }
  if (t == "Ntot") {
    int64_t i = op[0];
    return OpSum(Op("Nup", i)) + OpSum(Op("Ndn", i));
  }

  if (t == "NupNup") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nup", i) * Op("Nup", j));
  }
  if (t == "NdnNdn") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ndn", i) * Op("Ndn", j));
  }
  if (t == "NupNdn") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nup", i) * Op("Ndn", j));
  }
  if (t == "NdnNup") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ndn", i) * Op("Nup", j));
  }
  if (t == "NtotNtot") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += Op("Nup", i) * Op("Nup", j);
    r += Op("Nup", i) * Op("Ndn", j);
    r += Op("Ndn", i) * Op("Nup", j);
    r += Op("Ndn", i) * Op("Ndn", j);
    return r;
  }
  if (t == "SzSz") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.25 * (Op("Nup", i) * Op("Nup", j));
    r += -0.25 * (Op("Nup", i) * Op("Ndn", j));
    r += -0.25 * (Op("Ndn", i) * Op("Nup", j));
    r += 0.25 * (Op("Ndn", i) * Op("Ndn", j));
    return r;
  }
  if (t == "tJSzSz") {
    // tJSzSz = SzSz - 1/4 n n = -1/2 (Nup_i Ndn_j + Ndn_i Nup_j)
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -0.5 * (Op("Nup", i) * Op("Ndn", j));
    r += -0.5 * (Op("Ndn", i) * Op("Nup", j));
    return r;
  }

  if (t == "Sz") {
    int64_t i = op[0];
    return 0.5 * Op("Nup", i) + (-0.5) * Op("Ndn", i);
  }
  if (t == "S+") {
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cdn", i));
  }
  if (t == "S-") {
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cup", i));
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
  if (t == "SdotS") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("SzSz", {i, j})) + OpSum(Op("Exchange", {i, j}));
  }
  if (t == "tJSdotS") {
    int64_t i = op[0], j = op[1];
    return OpSum(Op("tJSzSz", {i, j})) + OpSum(Op("Exchange", {i, j}));
  }
  if (t == "Exchange") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("Cdagup", i) * Op("Cdn", i) * Op("Cdagdn", j) * Op("Cup", j));
    r += 0.5 * (Op("Cdagdn", i) * Op("Cup", i) * Op("Cdagup", j) * Op("Cdn", j));
    return r;
  }
  if (t == "ExchangeAsym") {
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("Cdagup", i) * Op("Cdn", i) * Op("Cdagdn", j) * Op("Cup", j));
    r -= 0.5 * (Op("Cdagdn", i) * Op("Cup", i) * Op("Cdagup", j) * Op("Cdn", j));
    return r;
  }

  if (t == "TotalN") {
    OpSum r;
    for (int64_t i = 0; i < alg.nsites; ++i) {
      r += Op("Ntot", i);
    }
    return r;
  }
  if (t == "TotalNup") {
    OpSum r;
    for (int64_t i = 0; i < alg.nsites; ++i) {
      r += Op("Nup", i);
    }
    return r;
  }
  if (t == "TotalNdn") {
    OpSum r;
    for (int64_t i = 0; i < alg.nsites; ++i) {
      r += Op("Ndn", i);
    }
    return r;
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
bool classify(std::string const &t, int &sector, int &creation) {
  if (t == "Cdagup") {
    sector = 0;
    creation = 0;
  } else if (t == "Cup") {
    sector = 0;
    creation = 1;
  } else if (t == "Cdagdn") {
    sector = 1;
    creation = 0;
  } else if (t == "Cdn") {
    sector = 1;
    creation = 1;
  } else {
    return false;
  }
  return true;
}
} // namespace

std::optional<OpSum> tj_simplify_symmetry(Monomial const &mono,
                                          Algebra const &alg) {
  (void)alg;
  int64_t n = mono.size();

  // Pass 1: same-site projected relations (site-major: dn before up).
  for (int64_t k = 0; k + 1 < n; ++k) {
    Op const &a = mono[k];
    Op const &b = mono[k + 1];
    if (a.size() != 1 || b.size() != 1 || a[0] != b[0]) {
      continue;
    }
    int sa, ca, sb, cb;
    if (!classify(a.type(), sa, ca) || !classify(b.type(), sb, cb)) {
      continue;
    }
    int64_t i = a[0];
    std::string ta = a.type(), tb = b.type();
    if (ca == cb && sa == sb) {
      return replace_pair(mono, k, OpSum{}); // nilpotent
    }
    if (ta == "Cup" && tb == "Cdagup") {
      // {Cup, Cdagup} = Id - Ndn
      return replace_pair(mono, k,
                          OpSum(Op("Id")) - OpSum(Op("Cdagup", i) * Op("Cup", i)) -
                              OpSum(Op("Cdagdn", i) * Op("Cdn", i)));
    }
    if (ta == "Cdn" && tb == "Cdagdn") {
      return replace_pair(mono, k,
                          OpSum(Op("Id")) - OpSum(Op("Cdagdn", i) * Op("Cdn", i)) -
                              OpSum(Op("Cdagup", i) * Op("Cup", i)));
    }
    if ((ta == "Cdn" && tb == "Cdagup") || (ta == "Cup" && tb == "Cdagdn")) {
      return replace_pair(mono, k, OpSum{}); // 0 (no double occupancy)
    }
    if (ta == "Cdagup" && tb == "Cdagdn") {
      return replace_pair(mono, k, -1.0 * (Op("Cdagdn", i) * Op("Cdagup", i)));
    }
    if (ta == "Cup" && tb == "Cdn") {
      return replace_pair(mono, k, -1.0 * (Op("Cdn", i) * Op("Cup", i)));
    }
  }

  // Pass 2: generic site-major sort.
  return sort_sites_step(mono, alg);
}

std::optional<OpSum> tj_simplify(Monomial const &mono, Algebra const &alg) {
  (void)alg;
  int64_t n = mono.size();
  for (int64_t k = 0; k + 1 < n; ++k) {
    Op const &a = mono[k];
    Op const &b = mono[k + 1];
    if (a.size() != 1 || b.size() != 1) {
      continue;
    }
    int sec_a, cre_a, sec_b, cre_b;
    if (!classify(a.type(), sec_a, cre_a) ||
        !classify(b.type(), sec_b, cre_b)) {
      continue;
    }
    int64_t sa = a[0], sb = b[0];

    if (sa != sb) {
      // Different sites: creation-major key (creation, sector, site).
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
        if (sec_a > sec_b) {
          return swap_pair(mono, k, -1.0); // order up(0) before dn(1)
        }
        continue;
      } else if ((cre_a == 0) && (cre_b == 1)) {
        continue; // Cdag then C, opposite spin: S+ / S- kept as-is
      } else {
        return replace_pair(mono, k, OpSum{}); // C then Cdag, opposite spin: 0
      }
    }
  }
  return std::nullopt;
}

} // namespace xdiag::algebra
