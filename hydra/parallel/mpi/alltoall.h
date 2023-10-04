#pragma once
#ifdef HYDRA_USE_MPI
#include <mpi.h>

namespace hydra::mpi {

template <class coeff_t>
int Alltoall(coeff_t *sendbuf, int sendcount, coeff_t *recvbuf, int recvcount,
             MPI_Comm comm);

template <class coeff_t>
int Alltoallv(coeff_t *sendbuf, int *sendcounts, int *sdispls, coeff_t *recvbuf,
              int *recvcounts, int *rdispls, MPI_Comm comm);

} // namespace hydra::mpi
#endif
