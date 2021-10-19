#pragma once

#include <mpi.h>
#include <hydra/mpi/datatype.h>

namespace hydra::mpi {

template <class TCoeffs>
int Allreduce(const TCoeffs *sendbuf, TCoeffs *recvbuf, int count, MPI_Op op,
              MPI_Comm comm);

} // namespace hydra::mpi
