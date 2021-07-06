#include "allreduce.h"

namespace hydra::mpi {

template <class TCoeffs>
int Allreduce(const TCoeffs *sendbuf, TCoeffs *recvbuf, int count, MPI_Op op,
              MPI_Comm comm) {
  MPI_Datatype type = datatype<TCoeffs>();
  return MPI_Allreduce(sendbuf, recvbuf, count, type, op, comm);
}

template int Allreduce<char>(const char *, char *, int, MPI_Op, MPI_Comm);
template int Allreduce<short>(const short *, short *, int, MPI_Op, MPI_Comm);
template int Allreduce<int>(const int *, int *, int, MPI_Op, MPI_Comm);
template int Allreduce<long>(const long *, long *, int, MPI_Op, MPI_Comm);
template int Allreduce<long long>(const long long *, long long *, int, MPI_Op,
                                  MPI_Comm);
template int Allreduce<unsigned char>(const unsigned char *, unsigned char *,
                                      int, MPI_Op, MPI_Comm);
template int Allreduce<unsigned short>(const unsigned short *, unsigned short *,
                                       int, MPI_Op, MPI_Comm);
template int Allreduce<unsigned int>(const unsigned int *, unsigned int *, int,
                                     MPI_Op, MPI_Comm);
template int Allreduce<unsigned long>(const unsigned long *, unsigned long *,
                                      int, MPI_Op, MPI_Comm);
template int Allreduce<unsigned long long>(const unsigned long long *,
                                           unsigned long long *, int, MPI_Op,
                                           MPI_Comm);
template int Allreduce<float>(const float *, float *, int, MPI_Op, MPI_Comm);
template int Allreduce<double>(const double *, double *, int, MPI_Op, MPI_Comm);

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
