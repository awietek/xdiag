#include "alltoall.h"

#include <hydra/mpi/datatype.h>

namespace hydra::mpi {

template <class TCoeffs>
int Alltoall(const TCoeffs *sendbuf, int sendcount, TCoeffs *recvbuf,
             int recvcount, MPI_Comm comm) {
  MPI_Datatype type = mpi::datatype<TCoeffs>();
  return MPI_Alltoall(sendbuf, sendcount, type, recvbuf, recvcount, type, comm);
}

template <class TCoeffs>
int Alltoallv(const TCoeffs *sendbuf, int *sendcounts, int *sdispls,
              TCoeffs *recvbuf, int *recvcounts, int *rdispls, MPI_Comm comm) {
  MPI_Datatype type = mpi::datatype<TCoeffs>();
  return MPI_Alltoallv(sendbuf, sendcounts, sdispls, type, recvbuf, recvcounts,
                       rdispls, type, comm);
}

// Special implementation for complex numbers
template <>
int Alltoallv<scomplex>(const scomplex *sendbuf, int *sendcounts, int *sdispls,
                        scomplex *recvbuf, int *recvcounts, int *rdispls,
                        MPI_Comm comm) {
  int n_mpi_tasks, ret;
  MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_tasks);
  for (int i = 0; i < n_mpi_tasks; ++i) {
    sendcounts[i] <<= 1;
    sdispls[i] <<= 1;
    recvcounts[i] <<= 1;
    rdispls[i] <<= 1;
  }
  ret = MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_FLOAT, recvbuf,
                      recvcounts, rdispls, MPI_FLOAT, comm);
  for (int i = 0; i < n_mpi_tasks; ++i) {
    sendcounts[i] >>= 1;
    sdispls[i] >>= 1;
    recvcounts[i] >>= 1;
    rdispls[i] >>= 1;
  }
  return ret;
}

template <>
int Alltoallv<complex>(const complex *sendbuf, int *sendcounts, int *sdispls,
                       complex *recvbuf, int *recvcounts, int *rdispls,
                       MPI_Comm comm) {
  int n_mpi_tasks, ret;
  MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_tasks);
  for (int i = 0; i < n_mpi_tasks; ++i) {
    sendcounts[i] <<= 1;
    sdispls[i] <<= 1;
    recvcounts[i] <<= 1;
    rdispls[i] <<= 1;
  }
  ret = MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE, recvbuf,
                      recvcounts, rdispls, MPI_DOUBLE, comm);
  for (int i = 0; i < n_mpi_tasks; ++i) {
    sendcounts[i] >>= 1;
    sdispls[i] >>= 1;
    recvcounts[i] >>= 1;
    rdispls[i] >>= 1;
  }
  return ret;
}

} // namespace hydra::mpi
