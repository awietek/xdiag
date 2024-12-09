#include "hc.hpp"

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
  } else if (type == "CDAGUP") {
    return Op("CUP", op.sites());
  } else if (type == "CUP") {
    return Op("CDAGUP", op.sites());
  } else if (type == "CDAGDN") {
    return Op("CDN", op.sites());
  } else if (type == "CDN") {
    return Op("CDAGDN", op.sites());
  } else { // default: the Op is hermitian

    if (op.hassites()) {
      if (op.hasmatrix()) {
        return Op(op.type(), op.sites(), op.matrix().hc());
      } else {
        return Op(op.type(), op.sites());
      }
    } else { // no sites defined
      if (op.hasmatrix()) {
        return Op(op.type(), op.matrix().hc());
      } else {
        return Op(op.type());
      }
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum hc(OpSum const &ops) try {
  OpSum ops_hc;
  for (auto [cpl, op] : ops.plain()) {
    ops_hc += conj(cpl.scalar()) * hc(op);
  }
  return ops_hc;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
} // namespace xdiag
