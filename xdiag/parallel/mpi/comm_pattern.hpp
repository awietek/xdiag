// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_USE_MPI
#include <xdiag/operators/op.hpp>
#include <xdiag/parallel/mpi/communicator.hpp>

namespace xdiag::mpi {

class CommPattern {
public:
  CommPattern() = default;
  bool contains(Op const &op) const;
  Communicator const &operator[](Op const &op) const;
  void append(Op const& op, Communicator const& comm);
private:
  std::vector<Op> ops_;
  std::vector<Communicator> comms_;
};

} // namespace xdiag::mpi
#endif
