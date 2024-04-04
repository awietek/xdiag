#include "datatype.h"

namespace xdiag::mpi {

template <> MPI_Datatype datatype<char>() { return MPI_CHAR; }
template <> MPI_Datatype datatype<short>() { return MPI_SHORT; }
template <> MPI_Datatype datatype<int>() { return MPI_INT; }
template <> MPI_Datatype datatype<long>() { return MPI_LONG_LONG_INT; }
template <> MPI_Datatype datatype<long long>() { return MPI_LONG_LONG_INT; }
template <> MPI_Datatype datatype<unsigned char>() { return MPI_UNSIGNED_CHAR; }
template <> MPI_Datatype datatype<unsigned short>() {
  return MPI_UNSIGNED_SHORT;
}
template <> MPI_Datatype datatype<unsigned int>() { return MPI_UNSIGNED; }
template <> MPI_Datatype datatype<unsigned long>() { return MPI_LONG_LONG_INT; }
template <> MPI_Datatype datatype<unsigned long long>() {
  return MPI_LONG_LONG_INT;
}
template <> MPI_Datatype datatype<float>() { return MPI_FLOAT; }
template <> MPI_Datatype datatype<double>() { return MPI_DOUBLE; }

} // namespace xdiag::mpi
