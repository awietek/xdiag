// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "allreduce.hpp"

namespace xdiag::mpi {

template <class coeff_t>
int Allreduce(coeff_t *sendbuf, coeff_t *recvbuf, int count, MPI_Op op,
              MPI_Comm comm) {
  MPI_Datatype type = datatype<coeff_t>();
  return MPI_Allreduce(sendbuf, recvbuf, count, type, op, comm);
}

template int Allreduce<char>(char *, char *, int, MPI_Op, MPI_Comm);
template int Allreduce<short>(short *, short *, int, MPI_Op, MPI_Comm);
template int Allreduce<int>(int *, int *, int, MPI_Op, MPI_Comm);
template int Allreduce<long>(long *, long *, int, MPI_Op, MPI_Comm);
template int Allreduce<long long>(long long *, long long *, int, MPI_Op,
                                  MPI_Comm);
template int Allreduce<unsigned char>(unsigned char *, unsigned char *, int,
                                      MPI_Op, MPI_Comm);
template int Allreduce<unsigned short>(unsigned short *, unsigned short *, int,
                                       MPI_Op, MPI_Comm);
template int Allreduce<unsigned int>(unsigned int *, unsigned int *, int,
                                     MPI_Op, MPI_Comm);
template int Allreduce<unsigned long>(unsigned long *, unsigned long *, int,
                                      MPI_Op, MPI_Comm);
template int Allreduce<unsigned long long>(unsigned long long *,
                                           unsigned long long *, int, MPI_Op,
                                           MPI_Comm);
template int Allreduce<double>(double *, double *, int, MPI_Op, MPI_Comm);

// Special implementation for complex numbers
template <>
int Allreduce<complex>(complex *sendbuf, complex *recvbuf, int count, MPI_Op op,
                       MPI_Comm comm) {
  return MPI_Allreduce(sendbuf, recvbuf, count << 1, MPI_DOUBLE, op, comm);
}

} // namespace xdiag::mpi
