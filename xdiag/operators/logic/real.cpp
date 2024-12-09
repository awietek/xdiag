#include "real.hpp"

#include <xdiag/operators/logic/types.hpp>
#include <xdiag/operators/logic/valid.hpp>

#include <string>
#include <vector>

namespace xdiag {

bool isreal(Op const &op) try {
  check_valid(op);

  std::string type = op.type();
  if (op.hasmatrix()) {
    return op.matrix().isreal();
  } else {
    return is_real_type(type);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool isreal(OpSum const &ops) try {
  for (auto [cpl, op] : ops.plain()) {
    if (!isreal(cpl.scalar()) || !isreal(op)) {
      return false;
    }
  }
  return true;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
