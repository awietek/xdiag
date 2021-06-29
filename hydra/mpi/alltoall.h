#pragma once

#include <mpi.h>

namespace hydra::mpi {

template <class TCoeffs>
int Alltoall(const TCoeffs *sendbuf, int sendcount, TCoeffs *recvbuf,
             int recvcount, MPI_Comm comm);

}
