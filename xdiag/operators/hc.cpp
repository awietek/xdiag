// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "hc.hpp"

#include <map>
#include <set>
#include <string>

#include <xdiag/math/matrix.hpp>
#include <xdiag/operators/types.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

// Two primitives, transpose and conj. Hermitian conjugation is
// hc = conj(transpose(.)).
//
// conj(type{sites}) = s * type{sites}, with s = -1 if the Fock-basis matrix of
// the type is purely imaginary (is_real = false in the registry), else s = +1.
//
// transpose(type{sites}) = s * partner{sites}. The partner is the type itself,
// unless it is listed in transpose_partner (the transpose of S+ is S-, the
// transpose of Cdag is C, ...). The sign is s = -1 if the matrix is
// antisymmetric, else s = +1.
//
// Every type not listed in transpose_partner or antisymmetric_types is symmetric
// and its own partner, so transpose and conj leave it unchanged. Unknown types
// are rejected by operators::check_valid.
//
// clang-format off
static const std::map<std::string, std::string> transpose_partner = {
  {"S+",      "S-"     }, {"S-",      "S+"     },
  {"A",       "Adag"   }, {"Adag",    "A"      },
  {"C",       "Cdag"   }, {"Cdag",    "C"      },
  {"Cdagup",  "Cup"    }, {"Cup",     "Cdagup" },
  {"Cdagdn",  "Cdn"    }, {"Cdn",     "Cdagdn" },
};

// Types with an antisymmetric Fock-basis matrix. This includes the real
// anti-hermitian ops (ExchangeAsym, HopAsym, ...) and the imaginary ones Sy and
// ScalarChirality, which are hermitian and antisymmetric at the same time.
static const std::set<std::string> antisymmetric_types = {
  "ExchangeAsym", "HopAsym", "HopupAsym", "HopdnAsym", "Sy", "ScalarChirality"};
// clang-format on

OpSum transpose(Op const &op) try {
  operators::check_valid(op);

  if (op.type() == "Matrix") {
    return OpSum(Op("Matrix", op.sites(), op.matrix().transpose()));
  } else {
    auto it = transpose_partner.find(op.type());
    std::string partner =
        (it != transpose_partner.end()) ? it->second : op.type();
    double phase = antisymmetric_types.count(op.type()) ? -1.0 : 1.0;
    if (op.hassites()) {
      return OpSum(phase, Op(partner, op.sites()));
    } else {
      return OpSum(phase, Op(partner));
    }
  }
}
XDIAG_CATCH

OpSum transpose(Monomial const &mono) try {
  // transpose(A_1 * ... * A_n) = transpose(A_n) * ... * transpose(A_1)
  OpSum mono_transpose(Monomial{}); // multiplicative identity (empty product)
  for (int64_t i = mono.size() - 1; i >= 0; --i) {
    mono_transpose = mono_transpose * transpose(mono[i]);
  }
  return mono_transpose;
}
XDIAG_CATCH

OpSum transpose(OpSum const &ops) try {
  OpSum ops_transpose;
  for (auto const &[coeff, mono] : ops.plain()) {
    // transpose(coeff * mono) = coeff * transpose(mono): no conjugation of coeff
    ops_transpose += coeff.scalar() * transpose(mono);
  }
  return ops_transpose;
}
XDIAG_CATCH

OpSum conj(Op const &op) try {
  operators::check_valid(op);

  if (op.type() == "Matrix") {
    return OpSum(Op("Matrix", op.sites(), op.matrix().conj()));
  } else {
    double phase = is_real_type(op.type()) ? 1.0 : -1.0;
    if (op.hassites()) {
      return OpSum(phase, Op(op.type(), op.sites()));
    } else {
      return OpSum(phase, Op(op.type()));
    }
  }
}
XDIAG_CATCH

OpSum conj(Monomial const &mono) try {
  // conj(A_1 * ... * A_n) = conj(A_1) * ... * conj(A_n)
  OpSum mono_conj(Monomial{}); // multiplicative identity (empty product = 1)
  for (int64_t i = 0; i < mono.size(); ++i) {
    mono_conj = mono_conj * conj(mono[i]);
  }
  return mono_conj;
}
XDIAG_CATCH

OpSum conj(OpSum const &ops) try {
  OpSum ops_conj;
  for (auto const &[coeff, mono] : ops.plain()) {
    // conj(coeff * mono) = conj(coeff) * conj(mono)
    ops_conj += conj(coeff.scalar()) * conj(mono);
  }
  return ops_conj;
}
XDIAG_CATCH

// hc = conj(transpose(.)). conj is multiplicative and transpose reverses
// products, so this also gives hc(A_1 * ... * A_n) = hc(A_n) * ... * hc(A_1) and
// hc(coeff * mono) = conj(coeff) * hc(mono).

OpSum hc(Op const &op) try {
  return conj(transpose(op));
}
XDIAG_CATCH

OpSum hc(Monomial const &mono) try {
  return conj(transpose(mono));
}
XDIAG_CATCH

OpSum hc(OpSum const &ops) try {
  return conj(transpose(ops));
}
XDIAG_CATCH

OpSum adjoint(Op const &op) try {
  return hc(op);
}
XDIAG_CATCH

OpSum adjoint(Monomial const &mono) try {
  return hc(mono);
}
XDIAG_CATCH

OpSum adjoint(OpSum const &ops) try {
  return hc(ops);
}
XDIAG_CATCH

} // namespace xdiag
