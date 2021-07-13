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
int Alltoallv(const TCoeffs *sendbuf, const int *sendcounts, const int *sdispls,
              TCoeffs *recvbuf, const int *recvcounts, const int *rdispls,
              MPI_Comm comm) {
  MPI_Datatype type = mpi::datatype<TCoeffs>();
  return MPI_Alltoallv(sendbuf, sendcounts, sdispls, type, recvbuf, recvcounts,
                       rdispls, type, comm);
}

template int Alltoallv<char>(const char *sendbuf, const int *sendcounts,
                             const int *sdispls, char *recvbuf,
                             const int *recvcounts, const int *rdispls,
                             MPI_Comm comm);
template int Alltoallv<short>(const short *sendbuf, const int *sendcounts,
                              const int *sdispls, short *recvbuf,
                              const int *recvcounts, const int *rdispls,
                              MPI_Comm comm);
template int Alltoallv<int>(const int *sendbuf, const int *sendcounts,
                            const int *sdispls, int *recvbuf,
                            const int *recvcounts, const int *rdispls,
                            MPI_Comm comm);
template int Alltoallv<long>(const long *sendbuf, const int *sendcounts,
                             const int *sdispls, long *recvbuf,
                             const int *recvcounts, const int *rdispls,
                             MPI_Comm comm);
template int Alltoallv<long long>(const long long *sendbuf,
                                  const int *sendcounts, const int *sdispls,
                                  long long *recvbuf, const int *recvcounts,
                                  const int *rdispls, MPI_Comm comm);
template int Alltoallv<unsigned char>(const unsigned char *sendbuf,
                                      const int *sendcounts, const int *sdispls,
                                      unsigned char *recvbuf,
                                      const int *recvcounts, const int *rdispls,
                                      MPI_Comm comm);
template int Alltoallv<unsigned short>(const unsigned short *sendbuf,
                                       const int *sendcounts,
                                       const int *sdispls,
                                       unsigned short *recvbuf,
                                       const int *recvcounts,
                                       const int *rdispls, MPI_Comm comm);
template int Alltoallv<unsigned int>(const unsigned int *sendbuf,
                                     const int *sendcounts, const int *sdispls,
                                     unsigned int *recvbuf,
                                     const int *recvcounts, const int *rdispls,
                                     MPI_Comm comm);
template int Alltoallv<unsigned long>(const unsigned long *sendbuf,
                                      const int *sendcounts, const int *sdispls,
                                      unsigned long *recvbuf,
                                      const int *recvcounts, const int *rdispls,
                                      MPI_Comm comm);
template int Alltoallv<unsigned long long>(const unsigned long long *sendbuf,
                                           const int *sendcounts,
                                           const int *sdispls,
                                           unsigned long long *recvbuf,
                                           const int *recvcounts,
                                           const int *rdispls, MPI_Comm comm);
template int Alltoallv<float>(const float *sendbuf, const int *sendcounts,
                              const int *sdispls, float *recvbuf,
                              const int *recvcounts, const int *rdispls,
                              MPI_Comm comm);
template int Alltoallv<double>(const double *sendbuf, const int *sendcounts,
                               const int *sdispls, double *recvbuf,
                               const int *recvcounts, const int *rdispls,
                               MPI_Comm comm);

// Special implementation for complex numbers
template <>
int Alltoallv<scomplex>(const scomplex *sendbuf, const int *sendcounts,
                        const int *sdispls, scomplex *recvbuf,
                        const int *recvcounts, const int *rdispls,
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
                       MPI_FLOAT, recvbuf, recvcounts_2.data(),
                       rdispls_2.data(), MPI_FLOAT, comm);
}

template <>
int Alltoallv<complex>(const complex *sendbuf, const int *sendcounts,
                       const int *sdispls, complex *recvbuf,
                       const int *recvcounts, const int *rdispls,
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
