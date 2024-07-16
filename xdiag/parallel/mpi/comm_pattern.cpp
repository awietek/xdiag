#include "comm_pattern.hpp"

namespace xdiag::mpi {

bool CommPattern::contains(Op const &op) {
  return comm_for_op_.count(op);
}

Communicator const &CommPattern::operator[](Op const &op) const try {
  if (contains(op)) {
    return comm_for_op_.at(op);
  } else {
    XDIAG_THROW("Cannot find communicator for Op");
    return Communicator();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

  Communicator &CommPattern::operator[](Op const &op try {
  if (contains(op)) {
    return comm_for_op_[op];
  } else {
    XDIAG_THROW("Cannot find communicator for op");
    return Communicator();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Communicator &operator[](Op const &op);

  } // namespace xdiag::mpi
