#include "alltoall.h"

#include <hydra/parallel/mpi/datatype.h>

namespace hydra::mpi {

///////////////////////////////////////////
// Alltoall
template <class coeff_t>
int Alltoall(coeff_t *sendbuf, int sendcount, coeff_t *recvbuf, int recvcount,
             MPI_Comm comm) {
  MPI_Datatype type = mpi::datatype<coeff_t>();
  return MPI_Alltoall(sendbuf, sendcount, type, recvbuf, recvcount, type, comm);
}

template int Alltoall<char>(char *sendbuf, int sendcount, char *recvbuf,
                            int recvcount, MPI_Comm comm);
template int Alltoall<short>(short *sendbuf, int sendcount, short *recvbuf,
                             int recvcount, MPI_Comm comm);
template int Alltoall<int>(int *sendbuf, int sendcount, int *recvbuf,
                           int recvcount, MPI_Comm comm);
template int Alltoall<long>(long *sendbuf, int sendcount, long *recvbuf,
                            int recvcount, MPI_Comm comm);
template int Alltoall<long long>(long long *sendbuf, int sendcount,
                                 long long *recvbuf, int recvcount,
                                 MPI_Comm comm);
template int Alltoall<unsigned char>(unsigned char *sendbuf, int sendcount,
                                     unsigned char *recvbuf, int recvcount,
                                     MPI_Comm comm);
template int Alltoall<unsigned short>(unsigned short *sendbuf, int sendcount,
                                      unsigned short *recvbuf, int recvcount,
                                      MPI_Comm comm);
template int Alltoall<unsigned int>(unsigned int *sendbuf, int sendcount,
                                    unsigned int *recvbuf, int recvcount,
                                    MPI_Comm comm);
template int Alltoall<unsigned long>(unsigned long *sendbuf, int sendcount,
                                     unsigned long *recvbuf, int recvcount,
                                     MPI_Comm comm);
template int Alltoall<unsigned long long>(unsigned long long *sendbuf,
                                          int sendcount,
                                          unsigned long long *recvbuf,
                                          int recvcount, MPI_Comm comm);

template int Alltoall<double>(double *sendbuf, int sendcount, double *recvbuf,
                              int recvcount, MPI_Comm comm);


template <>
int Alltoall<complex>(complex *sendbuf, int sendcounts, complex *recvbuf,
                      int recvcounts, MPI_Comm comm) {
  int n_mpi_tasks, ret;
  MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_tasks);
  sendcounts <<= 1;
  recvcounts <<= 1;
  ret = MPI_Alltoall(sendbuf, sendcounts, MPI_DOUBLE, recvbuf, recvcounts,
                     MPI_DOUBLE, comm);
  sendcounts >>= 1;
  recvcounts >>= 1;
  return ret;
}

///////////////////////////////////////////
// Alltoallv
template <class coeff_t>
int Alltoallv(coeff_t *sendbuf, int *sendcounts, int *sdispls, coeff_t *recvbuf,
              int *recvcounts, int *rdispls, MPI_Comm comm) {
  MPI_Datatype type = mpi::datatype<coeff_t>();
  return MPI_Alltoallv(sendbuf, sendcounts, sdispls, type, recvbuf, recvcounts,
                       rdispls, type, comm);
}

template int Alltoallv<char>(char *sendbuf, int *sendcounts, int *sdispls,
                             char *recvbuf, int *recvcounts, int *rdispls,
                             MPI_Comm comm);
template int Alltoallv<short>(short *sendbuf, int *sendcounts, int *sdispls,
                              short *recvbuf, int *recvcounts, int *rdispls,
                              MPI_Comm comm);
template int Alltoallv<int>(int *sendbuf, int *sendcounts, int *sdispls,
                            int *recvbuf, int *recvcounts, int *rdispls,
                            MPI_Comm comm);
template int Alltoallv<long>(long *sendbuf, int *sendcounts, int *sdispls,
                             long *recvbuf, int *recvcounts, int *rdispls,
                             MPI_Comm comm);
template int Alltoallv<long long>(long long *sendbuf, int *sendcounts,
                                  int *sdispls, long long *recvbuf,
                                  int *recvcounts, int *rdispls, MPI_Comm comm);
template int Alltoallv<unsigned char>(unsigned char *sendbuf, int *sendcounts,
                                      int *sdispls, unsigned char *recvbuf,
                                      int *recvcounts, int *rdispls,
                                      MPI_Comm comm);
template int Alltoallv<unsigned short>(unsigned short *sendbuf, int *sendcounts,
                                       int *sdispls, unsigned short *recvbuf,
                                       int *recvcounts, int *rdispls,
                                       MPI_Comm comm);
template int Alltoallv<unsigned int>(unsigned int *sendbuf, int *sendcounts,
                                     int *sdispls, unsigned int *recvbuf,
                                     int *recvcounts, int *rdispls,
                                     MPI_Comm comm);
template int Alltoallv<unsigned long>(unsigned long *sendbuf, int *sendcounts,
                                      int *sdispls, unsigned long *recvbuf,
                                      int *recvcounts, int *rdispls,
                                      MPI_Comm comm);
template int Alltoallv<unsigned long long>(unsigned long long *sendbuf,
                                           int *sendcounts, int *sdispls,
                                           unsigned long long *recvbuf,
                                           int *recvcounts, int *rdispls,
                                           MPI_Comm comm);

template int Alltoallv<double>(double *sendbuf, int *sendcounts, int *sdispls,
                               double *recvbuf, int *recvcounts, int *rdispls,
                               MPI_Comm comm);

// Special implementation for complex numbers
template <>
int Alltoallv<complex>(complex *sendbuf, int *sendcounts, int *sdispls,
                       complex *recvbuf, int *recvcounts, int *rdispls,
                       MPI_Comm comm) {
  int n_mpi_tasks;
  MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_tasks);
  std::vector<int> sendcounts_2(n_mpi_tasks, 0);
  std::vector<int> sdispls_2(n_mpi_tasks, 0);
  std::vector<int> recvcounts_2(n_mpi_tasks, 0);
  std::vector<int> rdispls_2(n_mpi_tasks, 0);
  for (int i = 0; i < n_mpi_tasks; ++i) {
    sendcounts_2[i] = sendcounts[i] * 2;
    sdispls_2[i] = sdispls[i] * 2;
    recvcounts_2[i] = recvcounts[i] * 2;
    rdispls_2[i] = rdispls[i] * 2;
  }
  return MPI_Alltoallv(sendbuf, sendcounts_2.data(), sdispls_2.data(),
                       MPI_DOUBLE, recvbuf, recvcounts_2.data(),
                       rdispls_2.data(), MPI_DOUBLE, comm);
}

} // namespace hydra::mpi
