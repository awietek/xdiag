#pragma once
#ifdef HYDRA_ENABLE_MPI

#include <mpi.h>
#include <hydra/mpi/datatype.h>

namespace hydra::mpi {

template <class TCoeffs>
int Allreduce(TCoeffs *sendbuf, TCoeffs *recvbuf, int count, MPI_Op op,
              MPI_Comm comm);

} // namespace hydra::mpi
#endif
