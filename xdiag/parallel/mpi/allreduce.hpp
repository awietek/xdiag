#pragma once
#ifdef XDIAG_USE_MPI

#include <mpi.h>

#include <xdiag/parallel/mpi/datatype.hpp>

namespace xdiag::mpi {

template <class coeff_t>
int Allreduce(coeff_t *sendbuf, coeff_t *recvbuf, int count, MPI_Op op,
              MPI_Comm comm);

} // namespace xdiag::mpi
#endif
