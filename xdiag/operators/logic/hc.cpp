// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "hc.hpp"

#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/types.hpp>
#include <xdiag/operators/logic/valid.hpp>

namespace xdiag {
Op hc(Op const &op) try {
  check_valid(op);

  std::string type = op.type();
  if (type == "S+") {
    return Op("S-", op.sites());
  } else if (type == "S-") {
    return Op("S+", op.sites());
  } else if (type == "Cdagup") {
    return Op("Cup", op.sites());
  } else if (type == "Cup") {
    return Op("Cdagup", op.sites());
  } else if (type == "Cdagdn") {
    return Op("Cdn", op.sites());
  } else if (type == "Cdn") {
    return Op("Cdagdn", op.sites());
  } else { // default: the type does not change

    if (op.hassites()) {
      if (op.hasmatrix()) {
        return Op(op.type(), op.sites(), op.matrix().hc());
      } else {
        return Op(op.type(), op.sites());
      }
    } else { // no sites defined
      return Op(op.type());
    }
  }
}
XDIAG_CATCH

OpSum hc(OpSum const &ops) try {
  OpSum ops_hc;
  for (auto [cpl, op] : ops.plain()) {
    std::string type = op.type();
    if ((type == "Exchange") || (type == "Hop") || (type == "Hopup") ||
        (type == "Hopdn")) {
      ops_hc += cpl.scalar() * hc(op);
    } else {
      ops_hc += conj(cpl.scalar()) * hc(op);
    }
  }
  return ops_hc;
}
XDIAG_CATCH

bool ishermitian(Op const &op) { return isapprox(op, hc(op)); }
bool ishermitian(OpSum const &ops) { return isapprox(ops, hc(ops)); }

} // namespace xdiag
