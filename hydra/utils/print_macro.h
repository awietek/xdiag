#pragma once

#include <hydra/utils/print.h>

#ifdef HYDRA_USE_MPI

#include <mpi.h>
namespace hydra::utils {
template <typename T> inline void print_pretty_mpi(const char *id, T X) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    hydra::utils::print_pretty(id, X);
  }
}

} // namespace hydra::utils
#define HydraPrint(X) hydra::utils::print_pretty_mpi(#X, X)

#else
#define HydraPrint(X) hydra::utils::print_pretty(#X, X)
#endif
