#include "alltoall.h"

#include <hydra/mpi/datatype.h>

namespace hydra::mpi {

///////////////////////////////////////////
// Alltoall
template <class TCoeffs>
int Alltoall(const TCoeffs *sendbuf, int sendcount, TCoeffs *recvbuf,
             int recvcount, MPI_Comm comm) {
  MPI_Datatype type = mpi::datatype<TCoeffs>();
  return MPI_Alltoall(sendbuf, sendcount, type, recvbuf, recvcount, type, comm);
}

template int Alltoall<char>(const char *sendbuf, int sendcount, char *recvbuf,
                            int recvcount, MPI_Comm comm);
template int Alltoall<short>(const short *sendbuf, int sendcount,
                             short *recvbuf, int recvcount, MPI_Comm comm);
template int Alltoall<int>(const int *sendbuf, int sendcount, int *recvbuf,
                           int recvcount, MPI_Comm comm);
template int Alltoall<long>(const long *sendbuf, int sendcount, long *recvbuf,
                            int recvcount, MPI_Comm comm);
template int Alltoall<long long>(const long long *sendbuf, int sendcount,
                                 long long *recvbuf, int recvcount,
                                 MPI_Comm comm);
template int Alltoall<unsigned char>(const unsigned char *sendbuf,
                                     int sendcount, unsigned char *recvbuf,
                                     int recvcount, MPI_Comm comm);
template int Alltoall<unsigned short>(const unsigned short *sendbuf,
                                      int sendcount, unsigned short *recvbuf,
                                      int recvcount, MPI_Comm comm);
template int Alltoall<unsigned int>(const unsigned int *sendbuf, int sendcount,
                                    unsigned int *recvbuf, int recvcount,
                                    MPI_Comm comm);
template int Alltoall<unsigned long>(const unsigned long *sendbuf,
                                     int sendcount, unsigned long *recvbuf,
                                     int recvcount, MPI_Comm comm);
template int Alltoall<unsigned long long>(const unsigned long long *sendbuf,
                                          int sendcount,
                                          unsigned long long *recvbuf,
                                          int recvcount, MPI_Comm comm);
template int Alltoall<float>(const float *sendbuf, int sendcount,
                             float *recvbuf, int recvcount, MPI_Comm comm);
template int Alltoall<double>(const double *sendbuf, int sendcount,
                              double *recvbuf, int recvcount, MPI_Comm comm);

// Special implementation for complex numbers
template <>
int Alltoall<scomplex>(const scomplex *sendbuf, int sendcounts,
                       scomplex *recvbuf, int recvcounts, MPI_Comm comm) {
  int n_mpi_tasks, ret;
  MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_tasks);
  sendcounts <<= 1;
  recvcounts <<= 1;
  ret = MPI_Alltoall(sendbuf, sendcounts, MPI_FLOAT, recvbuf, recvcounts,
                     MPI_FLOAT, comm);
  sendcounts >>= 1;
  recvcounts >>= 1;
  return ret;
}

template <>
int Alltoall<complex>(const complex *sendbuf, int sendcounts, complex *recvbuf,
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
template <class TCoeffs>
int Alltoallv(const TCoeffs *sendbuf, int *sendcounts, int *sdispls,
              TCoeffs *recvbuf, int *recvcounts, int *rdispls, MPI_Comm comm) {
  MPI_Datatype type = mpi::datatype<TCoeffs>();
  return MPI_Alltoallv(sendbuf, sendcounts, sdispls, type, recvbuf, recvcounts,
                       rdispls, type, comm);
}

template int Alltoallv<char>(const char *sendbuf, int *sendcounts, int *sdispls,
                             char *recvbuf, int *recvcounts, int *rdispls,
                             MPI_Comm comm);
template int Alltoallv<short>(const short *sendbuf, int *sendcounts,
                              int *sdispls, short *recvbuf, int *recvcounts,
                              int *rdispls, MPI_Comm comm);
template int Alltoallv<int>(const int *sendbuf, int *sendcounts, int *sdispls,
                            int *recvbuf, int *recvcounts, int *rdispls,
                            MPI_Comm comm);
template int Alltoallv<long>(const long *sendbuf, int *sendcounts, int *sdispls,
                             long *recvbuf, int *recvcounts, int *rdispls,
                             MPI_Comm comm);
template int Alltoallv<long long>(const long long *sendbuf, int *sendcounts,
                                  int *sdispls, long long *recvbuf,
                                  int *recvcounts, int *rdispls, MPI_Comm comm);
template int Alltoallv<unsigned char>(const unsigned char *sendbuf,
                                      int *sendcounts, int *sdispls,
                                      unsigned char *recvbuf, int *recvcounts,
                                      int *rdispls, MPI_Comm comm);
template int Alltoallv<unsigned short>(const unsigned short *sendbuf,
                                       int *sendcounts, int *sdispls,
                                       unsigned short *recvbuf, int *recvcounts,
                                       int *rdispls, MPI_Comm comm);
template int Alltoallv<unsigned int>(const unsigned int *sendbuf,
                                     int *sendcounts, int *sdispls,
                                     unsigned int *recvbuf, int *recvcounts,
                                     int *rdispls, MPI_Comm comm);
template int Alltoallv<unsigned long>(const unsigned long *sendbuf,
                                      int *sendcounts, int *sdispls,
                                      unsigned long *recvbuf, int *recvcounts,
                                      int *rdispls, MPI_Comm comm);
template int Alltoallv<unsigned long long>(const unsigned long long *sendbuf,
                                           int *sendcounts, int *sdispls,
                                           unsigned long long *recvbuf,
                                           int *recvcounts, int *rdispls,
                                           MPI_Comm comm);
template int Alltoallv<float>(const float *sendbuf, int *sendcounts,
                              int *sdispls, float *recvbuf, int *recvcounts,
                              int *rdispls, MPI_Comm comm);
template int Alltoallv<double>(const double *sendbuf, int *sendcounts,
                               int *sdispls, double *recvbuf, int *recvcounts,
                               int *rdispls, MPI_Comm comm);

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
