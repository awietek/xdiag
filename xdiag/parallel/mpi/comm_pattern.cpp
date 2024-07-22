#include "comm_pattern.hpp"

namespace xdiag::mpi {

bool CommPattern::contains(Op const &op) const {
  // Search is implemented pedestrian (instead of find) since no < is defined
  // for Op
  for (auto it = ops_.begin(); it != ops_.end(); ++it) {
    if (*it == op) {
      return true;
    }
  }
  return false;
}

Communicator const &CommPattern::operator[](Op const &op) const try {
  // Search is implemented pedestrian (instead of find) since no < is defined
  // for Op
  if (contains(op)) {
    int idx = 0;
    for (auto it = ops_.begin(); it != ops_.end(); ++it, ++idx) {
      if (*it == op) {
        return comms_[idx];
      }
    }
    XDIAG_THROW("Cannot find communicator for Op");
  } else {
    XDIAG_THROW("Cannot find communicator for Op");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void CommPattern::append(Op const &op, Communicator const &comm) {
  ops_.push_back(op);
  comms_.push_back(comm);
}

} // namespace xdiag::mpi
