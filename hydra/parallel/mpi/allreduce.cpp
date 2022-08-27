#include "allreduce.h"

namespace hydra::mpi {

template <class TCoeffs>
int Allreduce(TCoeffs *sendbuf, TCoeffs *recvbuf, int count, MPI_Op op,
              MPI_Comm comm) {
  MPI_Datatype type = datatype<TCoeffs>();
  return MPI_Allreduce(sendbuf, recvbuf, count, type, op, comm);
}

template int Allreduce<char>(char *, char *, int, MPI_Op, MPI_Comm);
template int Allreduce<short>(short *, short *, int, MPI_Op, MPI_Comm);
template int Allreduce<int>(int *, int *, int, MPI_Op, MPI_Comm);
template int Allreduce<long>(long *, long *, int, MPI_Op, MPI_Comm);
template int Allreduce<long long>(long long *, long long *, int, MPI_Op,
                                  MPI_Comm);
template int Allreduce<unsigned char>(unsigned char *, unsigned char *,
                                      int, MPI_Op, MPI_Comm);
template int Allreduce<unsigned short>(unsigned short *, unsigned short *,
                                       int, MPI_Op, MPI_Comm);
template int Allreduce<unsigned int>(unsigned int *, unsigned int *, int,
                                     MPI_Op, MPI_Comm);
template int Allreduce<unsigned long>(unsigned long *, unsigned long *,
                                      int, MPI_Op, MPI_Comm);
template int Allreduce<unsigned long long>(unsigned long long *,
                                           unsigned long long *, int, MPI_Op,
                                           MPI_Comm);
template int Allreduce<float>(float *, float *, int, MPI_Op, MPI_Comm);
template int Allreduce<double>(double *, double *, int, MPI_Op, MPI_Comm);

// Special implementation for complex numbers
template <>
int Allreduce<scomplex>(scomplex *sendbuf, scomplex *recvbuf, int count,
                        MPI_Op op, MPI_Comm comm) {
  return MPI_Allreduce(sendbuf, recvbuf, count << 1, MPI_FLOAT, op, comm);
}

template <>
int Allreduce<complex>(complex *sendbuf, complex *recvbuf, int count,
                       MPI_Op op, MPI_Comm comm) {
  return MPI_Allreduce(sendbuf, recvbuf, count << 1, MPI_DOUBLE, op, comm);
}

} // namespace hydra::mpi
