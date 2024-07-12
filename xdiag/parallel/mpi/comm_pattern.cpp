#include "comm_pattern.hpp"

namespace xdiag::mpi {

bool CommPattern::contains(Bond const &bond) {
  return comm_for_bond_.count(bond);
}

Communicator const &CommPattern::operator[](Bond const &bond) const try {
  if (contains(bond)) {
    return comm_for_bond_.at(bond);
  } else {
    XDIAG_THROW("Cannot find communicator for bond");
    return Communicator();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

  Communicator &CommPattern::operator[](Bond const &bond try {
  if (contains(bond)) {
    return comm_for_bond_[bond];
  } else {
    XDIAG_THROW("Cannot find communicator for bond");
    return Communicator();
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Communicator &operator[](Bond const &bond);

  } // namespace xdiag::mpi
