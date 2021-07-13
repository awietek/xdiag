#pragma once

#include <mpi.h>

namespace hydra::mpi {

template <class TCoeffs>
int Alltoall(const TCoeffs *sendbuf, int sendcount, TCoeffs *recvbuf,
             int recvcount, MPI_Comm comm);

template <class TCoeffs>
int Alltoallv(const TCoeffs *sendbuf, const int *sendcounts, const int *sdispls,
              TCoeffs *recvbuf, const int *recvcounts, const int *rdispls, MPI_Comm comm);

} // namespace hydra::mpi
