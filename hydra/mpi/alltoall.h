#pragma once

#include <mpi.h>

namespace hydra::mpi {

template <class TCoeffs>
int Alltoall(const TCoeffs *sendbuf, int sendcount, TCoeffs *recvbuf,
             int recvcount, MPI_Comm comm);

template <class TCoeffs>
int Alltoallv(const TCoeffs *sendbuf, int *sendcounts, int *sdispls,
              TCoeffs *recvbuf, int *recvcounts, int *rdispls, MPI_Comm comm);

} // namespace hydra::mpi
