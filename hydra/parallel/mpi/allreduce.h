#pragma once
#ifdef HYDRA_USE_MPI

#include <mpi.h>
#include <hydra/parallel/mpi/datatype.h>

namespace hydra::mpi {

template <class coeff_t>
int Allreduce(coeff_t *sendbuf, coeff_t *recvbuf, int count, MPI_Op op,
              MPI_Comm comm);

} // namespace hydra::mpi
#endif
