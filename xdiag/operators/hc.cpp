// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "hc.hpp"

#include <xdiag/operators/types.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

Op hc(Op const &op) try {
  operators::check_valid(op);

  auto const &info = info_of_type(op.type());

  // "Matrix" type: same hc_partner but matrix must be conjugate-transposed
  if (op.type() == "Matrix") {
    return Op("Matrix", op.sites(), op.matrix().hc());
  }

  // All other types: hc_partner encodes the target type; sites are unchanged
  if (op.hassites()) {
    return Op(info.hc_partner, op.sites());
  } else {
    return Op(info.hc_partner);
  }
}
XDIAG_CATCH

OpSum hc(OpSum const &ops) try {
  OpSum ops_hc;
  for (auto const &[coeff, mono] : ops.plain()) {
    ops_hc += conj(coeff.scalar()) * mono.hc();
  }
  return ops_hc;
}
XDIAG_CATCH

} // namespace xdiag
