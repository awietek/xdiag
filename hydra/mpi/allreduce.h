#pragma once

#include <mpi.h>
#include <hydra/mpi/datatype.h>

namespace hydra::mpi {

template <class TCoeffs>
int Allreduce(const TCoeffs *sendbuf, TCoeffs *recvbuf, int count, MPI_Op op,
              MPI_Comm comm) {
  MPI_Datatype type = datatype<TCoeffs>();
  return MPI_Allreduce(sendbuf, recvbuf, count, type, op, comm);
}

} // namespace hydra::mpi
