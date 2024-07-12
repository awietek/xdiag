#pragma once
#ifdef XDIAG_USE_MPI
#include <map>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/parallel/mpi/communicator.hpp>

namespace xdiag::mpi {

class CommPattern {
public:
  CommPattern = default();
  bool contains(Bond const &bond);
  Communicator const &operator[](Bond const &bond) const;
  Communicator &operator[](Bond const &bond);

private:
  std::map<Bond, Communicator> comm_for_bond_;
};

} // namespace xdiag::mpi
#endif
