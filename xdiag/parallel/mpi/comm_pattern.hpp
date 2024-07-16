#pragma once
#ifdef XDIAG_USE_MPI
#include <map>
#include <xdiag/operators/op.hpp>
#include <xdiag/parallel/mpi/communicator.hpp>

namespace xdiag::mpi {

class CommPattern {
public:
  CommPattern = default();
  bool contains(Op const &op);
  Communicator const &operator[](Op const &op) const;
  Communicator &operator[](Op const &op);

private:
  std::map<Op, Communicator> comm_for_op_;
};

} // namespace xdiag::mpi
#endif
