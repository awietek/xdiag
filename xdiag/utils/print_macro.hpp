#pragma once

#include <xdiag/utils/print.hpp>

#ifdef XDIAG_USE_MPI

#include <mpi.h>
namespace xdiag::utils {
template <typename T> inline void print_pretty_mpi(const char *id, T X) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    xdiag::utils::print_pretty(id, X);
  }
}

} // namespace xdiag::utils
#define XDIAG_PRINT(X) xdiag::utils::print_pretty_mpi(#X, X)

#else
#define XDIAG_PRINT(X) xdiag::utils::print_pretty(#X, X)
#endif
