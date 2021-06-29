#include "allreduce.h"

namespace hydra::mpi {

// Special implementation for complex numbers
template <>
int Allreduce<scomplex>(const scomplex *sendbuf, scomplex *recvbuf, int count,
                        MPI_Op op, MPI_Comm comm) {
  return MPI_Allreduce(sendbuf, recvbuf, count << 1, MPI_FLOAT, op, comm);
}

template <>
int Allreduce<complex>(const complex *sendbuf, complex *recvbuf, int count,
                       MPI_Op op, MPI_Comm comm) {
  return MPI_Allreduce(sendbuf, recvbuf, count << 1, MPI_DOUBLE, op, comm);
}

} // namespace hydra::mpi
